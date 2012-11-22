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

static t_class *pod_tilde_class;

typedef struct _pod_tilde
{
    t_object    x_obj;
    t_sample    x_f;
    t_int       flag;
    t_outlet*   bang;
    t_float     o_a1, o_a2, o_b0, o_b1, o_b2;
    t_float     m_a1, m_a2, m_b0, m_b1, m_b2;
    t_sample*   signal;
    t_int       window_size;
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


static t_int* pod_tilde_perform(t_int* w)
{
    t_pod_tilde *x = (t_pod_tilde *)(w[1]);     // x is the reference to the data struct
    t_sample  *in1 =    (t_sample *)(w[2]);     // in1 is an array of input samples
    int          n =           (int)(w[3]);     // n is the number of samples passed to this function
    
    int thresh_print = 0;
    
    for (int i = 0; i < n; i++)
    {
        // These two functions filter the signal, they're broken up because
        // it might make more sense to just use the middle ear filter
        x->signal[i] = pod_tilde_outer_filter(x, in1[i]);
        x->signal[i] = pod_tilde_middle_filter(x, x->signal[i]);
        
        // This is just the most basic threshold onset detector
        if (x->signal[i] > 0.9)
            thresh_print = 1;
    }
    
    if (thresh_print == 1)
    {
        x->flag = thresh_print;
        pod_tilde_print(x);
        outlet_bang(x->bang);
        x->flag = 0;
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
}

static void* pod_tilde_new(t_floatarg window_size)
{
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
    x->window_size = window_size;           // There should be a check to see if this is a power of two
    x->signal = (t_sample *)t_getbytes(x->window_size * sizeof(t_sample));
    
    
    post("pod~ v.0.1 by Gregoire Tronel, Jay Clark, and Scott McCoid");
    post("window size: %i", x->window_size);
    
    return (void *)x; 
}

void pod_tilde_setup(void)
{
    pod_tilde_class = class_new(gensym("pod~"), (t_newmethod)pod_tilde_new, (t_method)pod_tilde_free, sizeof(t_pod_tilde), CLASS_DEFAULT, A_DEFFLOAT, 0);
    
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
    
}
