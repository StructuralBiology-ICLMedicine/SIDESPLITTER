
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

// Filter map thread arguments structure
typedef struct{
  r_mrc      *in;
  double    *out;
  int32_t   size;
  int32_t   step;
  int32_t thread;
} map_arg;

// Mask making thread argument structure
typedef struct{
  r_mrc     *out;
  double     rad;
  double     cen;
  int32_t   size;
  int32_t size_2;
  int32_t   step;
  int32_t thread;
} make_mask_arg;

void make_mask_thread(make_mask_arg *arg);
// Make mask at diameter
// pthread function

void add_map_thread(map_arg *arg);
// Add MRC map in to out
// pthread function

void apply_mask_thread(map_arg *arg);
// Multiply out by in elementwise
// pthread function
