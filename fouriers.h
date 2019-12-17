
/*                                                                         
 * Copyright 14/08/2019 - Dr. Christopher H. S. Aylett                     
 *                                                                         
 * This program is free software; you can redistribute it and/or modify    
 * it under the terms of version 3 of the GNU General Public License as    
 * published by the Free Software Foundation.                              
 *                                                                         
 * This program is distributed in the hope that it will be useful,         
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           
 * GNU General Public License for more details - YOU HAVE BEEN WARNED!     
 *                                                                         
 * Program: SIDESPLITTER V0.1                                               
 *                                                                         
 * Authors: Chris Aylett                                                   
 *          Colin Palmer                                                   
 *                                                                         
 */

// Inclusions
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>
#include <pthread.h>
#include <unistd.h>
#include <complex.h>
#include <fftw3.h>

// Add FFT thread arguments structure
typedef struct {
  fftw_complex  *in;
  fftw_complex *out;
  int32_t      size;
  int32_t      step;
  int32_t    thread;
} add_fft_arg;

// Calc FSC thread arguments structure
typedef struct{
  fftw_complex   *half1;
  fftw_complex   *half2;
  long double numerator;
  long double denomin_1;
  long double denomin_2;
  int32_t          size;
  int32_t          step;
  int32_t        thread;
} calc_fsc_arg;

// Filter map thread arguments structure
typedef struct{
  fftw_complex  *in;
  fftw_complex *out;
  double      hires;
  double      lores;
  double        dim;
  int32_t full_size;
  int32_t      full;
  int32_t      size;
  int32_t      step;
  int32_t    thread;
} filter_arg;

// Spectrum thread arguments structure
typedef struct{
  fftw_complex *in1;
  fftw_complex *in2;
  double      *out1;
  double      *out2;
  double       *sum;
  double       *sub;
  int32_t        *n;
  double        dim;
  int32_t full_size;
  int32_t      full;
  int32_t      size;
  int32_t      step;
  int32_t    thread;
} spec_arg;

void add_fft_thread(add_fft_arg *arg);
// Add FFT in to out
// pthread function

void bandpass_filter_thread(filter_arg *arg);
// Apply bandpass to in and writes to out
// List node specifies resolutions
// pthread function

void lowpass_filter_thread(filter_arg *arg);
// Butterworth lowpass from in to out
// List node specifies resolution
// pthread function

void calc_fsc_thread(calc_fsc_arg *arg);
// Calculate FSC over map
// pthread function

void get_spec_thread(spec_arg *arg);
// Calculate spectrum over map
// pthread function

void apply_spec_thread(spec_arg *arg);
// Calculate spectrum over map
// pthread function
