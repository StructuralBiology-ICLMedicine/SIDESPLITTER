
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
 * Program: SIDESPLITTER V1.2                                               
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

// Noise and signal power thread arguments structure
typedef struct{
  r_mrc    *mask;
  double    *in1;
  double    *in2;
  double   noise;
  double   count;
  int32_t   size;
  int32_t   step;
  int32_t thread;
  long double sigma;
} max_arg;

// Probabilistic correction thread arguments structure
typedef struct{
  double    *in1;
  double    *in2;
  double   *out1;
  double   *out2;
  double   *ori1;
  double   *ori2;
  double   noise;
  double     rcv;
  int32_t   size;
  int32_t   step;
  int32_t thread;
} ass_vox_arg;

void calc_max_noise_thread(max_arg *arg);
// Calculate noise and signal power
// pthread function

void assign_voxels_thread(ass_vox_arg *arg);
// Correct according to probability
// pthread function

void taper_voxels_thread(ass_vox_arg *arg);
// Correct according to probability
// pthread function
