//
//  pod~.h
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

//----- Method declarations --------//

//Initialization
void pod_tilde_setup(void);
static void* pod_tilde_new(t_floatarg window_size, t_floatarg hop_size);
static void pod_tilde_create_window(t_pod_tilde* x);
static void new_bark_bands(t_pod_tilde* x);
static void create_filterbank(t_pod_tilde* x);

//Perform
static t_int* pod_tilde_perform(t_int* w);

//Ear Filters
    //Outer Ear
static t_float pod_tilde_outer_filter(t_pod_tilde* x, t_sample in);
static t_float pod_tilde_middle_filter(t_pod_tilde* x, t_sample in);

    //Inner Ear
static void multiply_filterbank(t_pod_tilde* x);
static void condense_analysis(t_pod_tilde* x);
static void multiply_loudness(t_pod_tilde* x);

//Peak Picking Helper Functions
static t_float accumulate_bin_differences(t_pod_tilde* x);
static void iterate_bark_bins(t_pod_tilde* x);

//Utilities
static int isPowerOfTwo(unsigned int x);
static float halfwave_rectify(float value);

//Memory managment
static void free_bark_bands(t_pod_tilde* x);
static void pod_tilde_free(t_pod_tilde* x);

//DSP
static void pod_tilde_dsp(t_pod_tilde* x, t_signal** sp);

//User Input
static void pod_tilde_set_window_type(t_pod_tilde* x, t_float number);
static void pod_tilde_set_debounce_threshold(t_pod_tilde* x, t_float number);
static void pod_tilde_set_upper_threshold(t_pod_tilde* x, t_float number);
static void pod_tilde_set_lower_threshold(t_pod_tilde* x, t_float number);





