//
//  pod~.c
//  pod
//
//
// Copyright (C) 2012 Scott McCoid, Gregoire Tronel, Jay Clark
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial
// portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "m_pd.h"
#include <math.h>
#include <stdio.h>
#define BLOCK_SIZE 64
#define FS 44100.0
#define PI 3.14159265359
#define TWO_PI (2 * PI)
#define NUM_BARKS 24

#define NUM_BARK_FILTER_BUFS 2

// define bark limits and centers
t_int bark_lim[25] =  { 20, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500 };
//t_int bark_ctr[24] = { 50, 150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500, 10500, 13500 };
t_int bark_ctr[26] = {0, 50, 150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500, 10500, 13500, 15500};
t_float band_weightings[24] = { 0.7762, 0.6854, 0.6647, 0.6373, 0.6255, 0.6170, 0.6139, 0.6107, 0.6127, 0.6329, 0.6380, 0.6430, 0.6151, 0.6033, 0.5914, 0.5843, 0.5895, 0.5947, 0.6237, 0.6703, 0.6920, 0.7137, 0.7217, 0.7217 };


static t_class  *pod_tilde_class;

typedef struct _bark_bin
{
    t_float*    band;
    
} t_bark_bin;

typedef struct _pod_tilde
{
    t_object    x_obj;
    t_sample    x_f;
    t_outlet*   bang;
    t_outlet*   mag_outlet;
    t_outlet*   bin_diffs;
    t_float     o_a1, o_a2, o_b0, o_b1, o_b2;
    t_float     m_a1, m_a2, m_b0, m_b1, m_b2;
    t_sample*   signal;                         // this holds samples
    t_sample*   analysis;                       // this holds analysis values
    t_int       window_size;
    t_float*    window;                         
    t_int       window_type;
    t_int       hop_size;
    t_int       dsp_tick;
    t_int       half_window_size;
    
    //peak picking
    t_float     bark_bins[24];
    t_float     prev_bark_bins[24];
    t_float     u_threshold, l_threshold;
    t_float     bark_difference;
    t_float     peak_value;
    t_int       flag;
    t_int       debounce_iterator;
    t_int       debounce_threshold;
    
    // filterbank
    t_bark_bin  filter_bands[2];
    t_sample*   filtered_odd;
    t_sample*   filtered_even;
    
    
} t_pod_tilde;


static void pod_tilde_print(t_pod_tilde* x)
{
    post("Threshold: %i", (int)x->flag);
}

static void linspace(t_int low, t_int high, t_int subdiv, t_float* line_buffer)
{
    t_float iter = (high - low) / subdiv;
    line_buffer[0] = low;
    
    for (int i = 1; i < subdiv; i++)
        line_buffer[i] = line_buffer[i - 1] + iter;
}

static void multiply_filterbank(t_pod_tilde* x)
{
    for (int i = 0; i < x->half_window_size; i++)
    {
        x->filtered_odd[i] = x->filter_bands[0].band[i] * x->analysis[i];           // non overlapping bands starting at 0
        x->filtered_even[i] = x->filter_bands[1].band[i] * x->analysis[i];          // non overlapping bands starting at 50 (first bark center)
        x->analysis[i] = x->filtered_odd[i] + x->filtered_even[i];
    }
}

static void condense_analysis(t_pod_tilde* x)
{
    float period = FS / x->window_size;
    
    for (int i = 0; i < NUM_BARKS; i++)
    {
        for (int j = 0; j < x->half_window_size; j++)
        {
            float frequency = period * j;
            
            if (frequency >= bark_ctr[i] && frequency < bark_ctr[i + 2])
                x->bark_bins[i] += x->analysis[j] * band_weightings[i];
        }
    }
}


static void create_filterbank(t_pod_tilde* x)
{
    float period = FS / x->window_size;
    int direction = 0;                     // direction is either +1 for increasing or -1 for decreasing
    
    float length, slope, point;
    
    // NUM_BARKS is still 24, but we have an array of length 26, so we've added lower and upper limits
    for (int i = 0; i < NUM_BARKS; i++)
    {
        for (int j = 0; j < x->half_window_size; j++)
        {
            float frequency = period * j;
            
            if (frequency >= bark_ctr[i] && frequency < bark_ctr[i + 1])
            {
                direction = 1;
                length = bark_ctr[i + 1] - bark_ctr[i];
            }
            else if (frequency >= bark_ctr[i + 1] && frequency < bark_ctr[i + 2])
            {
                direction = -1;
                length = bark_ctr[i + 2] - bark_ctr[i + 1];
            }
            else
                direction = 0;  // this means we're over the bounds and don't want to deal with it

            if (direction != 0)
            {
                slope = direction / length;
                point = 1 - slope * bark_ctr[i + 1];
            
                x->filter_bands[i % 2].band[j] = slope * frequency + point;         // y = mx + b
            }
        }
    }
}

static t_float pod_tilde_outer_filter(t_pod_tilde* x, t_sample in)
{
    //a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb) - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
    static t_float o_x1 = 0.0;
    static t_float o_x2 = 0.0;
    static t_float o_y1 = 0.0;
    static t_float o_y2 = 0.0;
    t_float out = 0.0;
    
    out = x->o_b0 * in + x->o_b1 * o_x1 + x->o_b2 * o_x2 - x->o_a1 * o_y1 - x->o_a2 * o_y2;
    
    o_x2 = o_x1;
    o_x1 = in;
    o_y2 = o_y1;
    o_y1 = out;
    
    return out;
}

static t_float pod_tilde_middle_filter(t_pod_tilde* x, t_sample in)
{
    static t_float m_x1 = 0.0;
    static t_float m_x2 = 0.0;
    static t_float m_y1 = 0.0;
    static t_float m_y2 = 0.0;
    t_float out = 0.0;
    
    m_y1 = x->m_b0 * in + x->m_b1 * m_x1 + x->m_b2 * m_x2  - x->m_a2 * m_y2;
    out = x->m_a2 * m_y1;
    
    m_x2 = m_x1;
    m_x1 = in;
    m_y2 = m_y1;
    m_y1 = out;
        
    return out;
}

static int isPowerOfTwo(unsigned int x){
    //Complement and Compare
    return ((x != 0) && ((x & (~x + 1)) == x));
}

static void pod_tilde_set_window_type(t_pod_tilde* x, t_float number){
    
    int selection = (int) number;
    if (selection>=0 & selection<=1) {
        x->window_type = selection;
    }
    else (post("Invalid windowing parameter"));
    
}

static void pod_tilde_set_debounce_threshold(t_pod_tilde* x, t_float number){
    
    int selection = (int) number;
    // need to add error checking
    x->debounce_threshold = selection;
    
}

static void pod_tilde_set_upper_threshold(t_pod_tilde* x, t_float number){
    
    // need to add error checking
    x->u_threshold = number;
    
}

static void pod_tilde_set_lower_threshold(t_pod_tilde* x, t_float number){
    
    // need to add error checking
    x->l_threshold = number;
    
}

static t_float accumulate_bin_differences(t_pod_tilde* x){
    
    t_float diff = 0;
    t_int length = sizeof(x->bark_bins) / sizeof(t_float);
    for (int i = 0; i < length; i++){
        diff += x->bark_bins[i] - x->prev_bark_bins[i];
    }
    
    outlet_float(x->bin_diffs, diff);
    
    return diff;
}

static void iterate_bark_bins(t_pod_tilde* x){
    
    t_int length = sizeof(x->bark_bins) / sizeof(t_float);
    for (int i = 0; i < length; i++) {
        x->prev_bark_bins[i] = x->bark_bins[i];
    }
    
}


static t_int* pod_tilde_perform(t_int* w)
{
    t_pod_tilde *x = (t_pod_tilde *)(w[1]);     // x is the reference to the data struct
    t_sample  *in1 =    (t_sample *)(w[2]);     // in1 is an array of input samples
    int          n =           (int)(w[3]);     // n is the number of samples passed to this function
    
    int size_diff = x->window_size - n;
    
    // This takes part of the signal buffer and shifts it to the front
    for (int i = 0; i < size_diff; i++)
        x->signal[i] = x->signal[i + n];
    
    // This takes the new samples filters them and puts them into the buffer
    for (int i = 0; i < n; i++)
    {
        x->signal[size_diff + i] = pod_tilde_outer_filter(x, in1[i]);
        x->signal[size_diff + i] = pod_tilde_middle_filter(x, x->signal[size_diff + i]);
    }
   
    // Increase the dsp_tick variable by the number of samples passed to callback
    x->dsp_tick += n;
    
    // If the dsp_tick reaches the hop_size value, then we do our processing
    if (x->dsp_tick >= x->hop_size)
    {
        x->dsp_tick = 0;
        
        // do windowing
        for (int i = 0; i < x->window_size; i++)
            x->analysis[i] = x->signal[i] * x->window[i];           // analysis is windowed signal
        
        // take fft
        mayer_realfft(x->window_size, x->analysis);
        
        // Zero out frequencies at DC and Nyquist
        x->analysis[0] = 0.0;
        x->analysis[x->window_size / 2] = 0.0;
        
        // Get the magnitude and assign it to the first half of the analysis buffer
        for (int i = 0; i < x->window_size / 2; i++)
        {
            int i_index = x->window_size - i;
            x->analysis[i] = sqrt((x->analysis[i] * x->analysis[i]) + (x->analysis[i_index] * x->analysis[i_index]));
        }
    
        // multiply analysis buffer by the filterbank
        multiply_filterbank(x);
        
        // turn analysis buffer into summed bark bins (i.e. turn 1 x half_windowsize vector -> 1 x 24 vector)
        condense_analysis(x);
        
        
        // -- spectral flux peak picking -- //
        
        //check for initial case
        if (x->prev_bark_bins != NULL) {
        
        //subtract this frame from last to get to our feature space
        x->bark_difference = accumulate_bin_differences(x);

        
        //Is our flag raised?
        switch (x->flag) {
                
            case 0: //Flag is down.
                
                //Lets check if we're above the upper threshold
                if (x->bark_difference > x->u_threshold) {
                    
                    //Let's flag this spot for a potential onset and hang on to that peak value if it ends up being one
                    x->flag = 1;
                    x->debounce_iterator = 1;
                    x->peak_value=x->bark_difference;
                }

                //otherwise, we'll keep waiting for an onset. 
                
                break;
                
            case 1: //Flag is up.
                
                //did we go even higher above the threshold?
                if (x->bark_difference > x->peak_value) {
                    
                    //flag this as a better estimate for the onset.
                    x->flag = 1;
                    x->debounce_threshold = 1;
                    x->peak_value = x->bark_difference;
                    
                }
                
                //if not...
                else{
                    
                    //Have we gone beyond our debouncing window?
                    if (x->debounce_iterator > x->debounce_threshold) {
                        
                        //onset verified!
                        outlet_bang(x->bang);
                        outlet_float(x->mag_outlet, x->peak_value); 
            
                        x->debounce_iterator = 0;
                        x->flag = 0;
                        
                    }
                    
                    else{
                        
                        //are we below our lower threshold?
                        if(x->bark_difference < x->l_threshold){
                            
                            //onset verified!
                            outlet_bang(x->bang);
                            outlet_float(x->mag_outlet, x->peak_value);
                            
                            x->debounce_iterator = 0;
                            x->flag = 0;
                            
                            
                        }
                        
                        //we have a peak flagged, but we haven't increased or crossed the lower threshold yet.
                        //Lets wait a bit longer to make sure our tagged peak is an onset
                        else x->debounce_iterator++;
                    }
                      
                    
                }

                break;
        }
        }
        
        iterate_bark_bins(x);
        
    }
        
    return (w + 4);
}

static void pod_tilde_dsp(t_pod_tilde* x, t_signal** sp)
{
    dsp_add(pod_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

static void new_bark_bands(t_pod_tilde* x)
{
    for (int i = 0; i < NUM_BARK_FILTER_BUFS; i++)
    {
        x->filter_bands[i].band = (t_float *)t_getbytes((x->window_size / 2) * sizeof(t_float));
    
        for (int j = 0; j < x->half_window_size; j++)
            x->filter_bands[i].band[j] = 0.0;
    }
    
}

static void free_bark_bands(t_pod_tilde* x)
{
    for (int i = 0; i < NUM_BARK_FILTER_BUFS; i++)
        t_freebytes(x->filter_bands[i].band, (x->half_window_size) * sizeof(t_float));
}

static void pod_tilde_free(t_pod_tilde* x)
{
    t_freebytes(x->signal, x->window_size * sizeof(t_sample));
    t_freebytes(x->analysis, x->window_size * sizeof(t_sample));
    t_freebytes(x->window, x->window_size * sizeof(t_float));
    
    t_freebytes(x->filtered_odd, x->half_window_size * sizeof(t_sample));
    t_freebytes(x->filtered_even, x->half_window_size * sizeof(t_sample));
    
    free_bark_bands(x);
}

static void pod_tilde_create_window(t_pod_tilde* x)
{
    
    switch (x->window_type) {
        case 0:
            // Hanning
            for (int i = 0; i < x->window_size; i++)
                x->window[i] = 0.5 * (1 - cos((TWO_PI * i) / (x->window_size - 1)));
            break;
            
        case 1:
            // Hamming
            for (int i = 0; i < x->window_size; i++)
                x->window[i] = 0.54 - 0.46 * (cos((TWO_PI * i) / (x->window_size - 1)));
        default:
            post("Unexpected windowing method");
            break;
    }
}

static void* pod_tilde_new(t_floatarg window_size, t_floatarg hop_size)
{
    
    post("pod~ v.0.1 by Gregoire Tronel, Jay Clark, and Scott McCoid");
    
    t_pod_tilde *x = (t_pod_tilde *)pd_new(pod_tilde_class);
        
    // Leftmost outlet outputs a bang
    x->bang = outlet_new(&x->x_obj, &s_bang);
    x->mag_outlet = outlet_new(&x->x_obj, &s_float);
    x->bin_diffs = outlet_new(&x->x_obj, &s_float);
    
    // Initialize filter coeffs
    // Outer
    x->o_a1 = 0.0;
    x->o_a2 = 0.0;
    x->o_b0 = 0.0;
    x->o_b1 = 0.7221;
    x->o_b2 = -0.6918;
    
    // Middle
    x->m_a1 = 1.6456;
    x->m_a2 = 0.6791;
    x->m_b0 = 0.8383;
    x->m_b1 = 0.0;
    x->m_b2 = -0.8383;
    
    // Window Size
    if (! isPowerOfTwo(window_size)){
        post("Window size must be a power of two. Applying default window.");
        x->window_size = 1024;
    }
    else x->window_size = window_size;
    
    x->half_window_size = x->window_size / 2;

    x->signal = (t_sample *)t_getbytes(x->window_size * sizeof(t_sample));
    x->analysis = (t_sample *)t_getbytes(x->window_size * sizeof(t_sample));
    x->window = (t_float *)t_getbytes(x->window_size * sizeof(t_float));
    x->filtered_odd = (t_sample *)t_getbytes(x->half_window_size * sizeof(t_sample));
    x->filtered_even= (t_sample *)t_getbytes(x->half_window_size * sizeof(t_sample));
    
        
    for (int i = 0; i < x->window_size; i++)
    {
        x->signal[i] = 0.0;
        x->analysis[i] = 0.0;
        x->window[i] = 0.0;
    }
    x->window_type = 0; //Default hanning
    pod_tilde_create_window(x);
    
    // create arrays for each bark filter band
    new_bark_bands(x);
    
    // create filter-bank associated with window size
    create_filterbank(x);
    
    if (! isPowerOfTwo(hop_size)){
        post("Hop size must be a power of two. Applying default hop size.");
        x->hop_size = 256;
    }
    else x->hop_size = hop_size; // This is in samples
    x->dsp_tick = 0;
    
    post("window size: %i", x->window_size);
    post("hop size: %i", x->hop_size);
    
    
    //Peak picking.
    //Lower our flag. No onsets yet!
    
    x->flag = 0;
    x->debounce_iterator=0;
    x->debounce_threshold=5;
    x->u_threshold = 1000;
    x->l_threshold = 900;

    
    return (void *)x; 
}





void pod_tilde_setup(void)
{
    pod_tilde_class = class_new(gensym("pod~"), (t_newmethod)pod_tilde_new, (t_method)pod_tilde_free, sizeof(t_pod_tilde), CLASS_DEFAULT, A_DEFFLOAT, A_DEFFLOAT, 0);
    
    CLASS_MAINSIGNALIN(pod_tilde_class, t_pod_tilde, x_f);
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_dsp,
        gensym("dsp"),
        0
    );
    
    class_addmethod(pod_tilde_class, (t_method)pod_tilde_print, gensym("print"), 0);
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_set_window_type,
        gensym("window"),
        A_FLOAT,
        0
    );
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_set_debounce_threshold,
        gensym("debounce"),
        A_FLOAT,
        0
    );
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_set_upper_threshold,
        gensym("upper"),
        A_FLOAT,
        0
        );
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_set_lower_threshold,
        gensym("lower"),
        A_FLOAT,
        0
        );
    

}
