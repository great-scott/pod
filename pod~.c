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
#define BLOCK_SIZE 64
#define PI 3.14159265359
#define TWO_PI (2 * PI)

static t_class *pod_tilde_class;

typedef struct _pod_tilde
{
    t_object    x_obj;
    t_sample    x_f;
    t_outlet*   bang;
    t_float     o_a1, o_a2, o_b0, o_b1, o_b2;
    t_float     m_a1, m_a2, m_b0, m_b1, m_b2;
    t_sample*   signal;                         // this holds samples
    t_sample*   analysis;                       // this holds analysis values
    t_int       window_size;
    t_float*    window;                         
    t_int       window_type;
    t_int       hop_size;
    t_int       dsp_tick;
    
    //onset detection
    t_float     bark_bins[24];
    t_float     prev_bark_bins[24];
    t_float     u_threshold, l_threshold;
    t_float     bark_difference;
    t_float     peak_value;
    t_int       flag;
    t_int       debounce_iterator;
    t_int       debounce_threshold;
    
    
    
} t_pod_tilde;


static void pod_tilde_print(t_pod_tilde* x)
{
    post("Threshold: %i", (int)x->flag);
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

static t_float accumulate_bin_differences(t_pod_tilde* x){
    
    t_float diff = 0;
    for (int i = 0; i < sizeof(x->bark_bins); i++){
        diff = diff + (x->bark_bins[i] - x->prev_bark_bins[i]);
    }
    
    return diff;
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
            x->analysis[i] = x->signal[i] * x->window[i];
        
        // take fft
        mayer_realfft(x->window_size, x->analysis);
        
        // Get the magnitude and assign it to the first half of the analysis buffer
        for (int i = 0; i < x->window_size / 2; i++)
        {
            int i_index = x->window_size - i;
            x->analysis[i] = sqrt((x->analysis[i] * x->analysis[i]) + (x->analysis[i_index] * x->analysis[i_index]));
        }
        
        // multiply analysis buffer by the filterbank
          
        
        
        
////////////////////////////// where the magic happens /////////////////////////////////////////////////
        
        
        
            //assume filterbanks are generalized into each bin of the length 24 array "bark_bins"
        
        //NEED TO ADDRESS INITIAL CASE
                
        //subtract this window from last to get to our feature space
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
                
                //can we go even higher above the threshold?
                if (x->bark_difference > x->peak_value) {
                    
                    //flag this as a better estimate for the onset.
                    x->flag = 1;
                    x->debounce_threshold = 1;
                    x->peak_value = x->bark_difference;
                    
                }
                
                else{
                    
                    //Have we gone beyond our debouncing window?
                    if (x->debounce_iterator > x->debounce_threshold) {
                        
                        //onset verified!
                        outlet_bang(x->bang);
                        //outlet_float(x, x->peak_value); //fix this after pull.
                        
                        x->debounce_iterator = 0;
                        x->flag = 0;
                        
                    }
                    
                    else{
                        
                        //are we below our lower threshold?
                        if(x->bark_difference < x->l_threshold){
                            
                            //onset verified!
                            outlet_bang(x->bang);
                            //outlet_float(x, x->peak_value); //fix this after pull.
                            
                            x->debounce_iterator = 0;
                            x->flag = 0;
                            
                            
                        }
                        
                        //we have a onset flagged, but we haven't increased or crossed the lower threshold yet.
                        //Lets wait a bit longer to make sure our tagged onset is legit
                        else x->debounce_iterator++;
                    
                    }
                    
                    
                    
                    
                }
                
                
              

                break;
        }
        
        
        
        
    }
        
    return (w + 4);
}

static void pod_tilde_dsp(t_pod_tilde* x, t_signal** sp)
{
    dsp_add(pod_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

static void pod_tilde_free(t_pod_tilde* x)
{
    t_freebytes(x->signal, x->window_size * sizeof(t_sample));
    t_freebytes(x->analysis, x->window_size * sizeof(t_sample));
    t_freebytes(x->window, x->window_size * sizeof(t_float));
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
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
        
    // Leftmost outlet outputs a bang
    x->bang = outlet_new(&x->x_obj, &s_bang);
    
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
        post("Window size must be a power of two. Applying defult window.");
        x->window_size = 1024;
    }
    else x->window_size = window_size;

    x->signal = (t_sample *)t_getbytes(x->window_size * sizeof(t_sample));
    x->analysis = (t_sample *)t_getbytes(x->window_size * sizeof(t_sample));
    x->window = (t_float *)t_getbytes(x->window_size * sizeof(t_float));
    for (int i = 0; i < x->window_size; i++)
    {
        x->signal[i] = 0.0;
        x->analysis[i] = 0.0;
        x->window[i] = 0.0;
    }
    x->window_type = 0; //Default hanning
    pod_tilde_create_window(x);
    
    if (! isPowerOfTwo(hop_size)){
        post("Hop size must be a power of two. Applying defult hop size.");
        x->hop_size = 256;
    }
    else x->hop_size = hop_size; // This is in samples
    x->dsp_tick = 0;
    
    post("window size: %i", x->window_size);
    post("hop size: %i", x->hop_size);
    
    
    
    //Lower our flag. No onsets yet!
    x->flag = 0;
    x->debounce_iterator=0;
    x->debounce_threshold=5;
    
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
        (t_method)pod_tilde_outer_filter,
        gensym("outer_filter"),
        A_DEFFLOAT,
        0
    );
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_middle_filter,
        gensym("middle_filter"),
        A_DEFFLOAT,
        0
    );
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_create_window,
        gensym("create_window"),
        0
    );
    
    class_addmethod(
        pod_tilde_class,
        (t_method)pod_tilde_set_window_type,
        gensym("window"),
        A_FLOAT,
        0
    );
    
}
