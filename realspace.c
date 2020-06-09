
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

// Library header inclusion for linking                                  
#include "sidesplitter.h"
#include "realspace.h"

// Make mask from radius in voxels
r_mrc *make_msk(r_mrc *in, double rad, int32_t nthreads){
  r_mrc *out = malloc(sizeof(r_mrc));
  int32_t size = in->n_crs[0], i;
  double cen = (double) size / 2;
  memcpy(out, in, sizeof(r_mrc));
  out->data = calloc(size * size * size, sizeof(float));
  pthread_t threads[nthreads];
  make_mask_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].rad = rad * rad;
    arg[i].out = out;
    arg[i].size = size;
    arg[i].size_2 = size * size;
    arg[i].cen = cen;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) make_mask_thread, &arg[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
  }
  return out;
}

void make_mask_thread(make_mask_arg *arg){
  double i, j, k, norm;
  for(int32_t _k = 0; _k < arg->size; _k++){
    k = (double) _k - arg->cen;
    k = k * k;
    for(int32_t _j = 0; _j < arg->size; _j++){
      j = (double) _j - arg->cen;
      j = j * j;
      for(int32_t _i = arg->thread; _i < arg->size; _i += arg->step){
        i = (double) _i - arg->cen;
        i = i * i;
        norm = (double) k + j + i;
        arg->out->data[ _k * arg->size_2 + _j * arg->size + _i ] = 1.0 / sqrt(1.0 + pow((norm / arg->rad), 8.0));
      }
    }
  }
  return;
}

// Add MRC map in to out
void add_map(r_mrc *in, double *out, int32_t nthreads){
  int32_t size = in->n_crs[0] * in->n_crs[1] * in->n_crs[2], i;
  pthread_t threads[nthreads];
  map_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) add_map_thread, &arg[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
  }
  return;
}

void add_map_thread(map_arg *arg){
  int32_t i;
  for (i = arg->thread; i < arg->size; i += arg->step){
    arg->out[i] += (double) arg->in->data[i];
  }
  return;
}

// Multiply out by in elementwise
void apply_mask(r_mrc *in, double *out, int32_t nthreads){
  int32_t size = in->n_crs[0] * in->n_crs[1] * in->n_crs[2], i;
  pthread_t threads[nthreads];
  map_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) apply_mask_thread, &arg[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
  }
  return;
}

void apply_mask_thread(map_arg *arg){
  int32_t i;
  for (i = arg->thread; i < arg->size; i += arg->step){
    arg->out[i] *= (double) arg->in->data[i];
  }
  return;
}
