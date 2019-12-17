
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

// Library header inclusion for linking                                  
#include "sidesplitter.h"
#include "suppress.h"

// Suppress noise between in/out
double suppress_noise(double *in1, double *in2, double *out1, double *out2, r_mrc *mask, list *node, int32_t size, int32_t nthreads){
  int32_t i, max = size * size * size;
  pthread_t threads[nthreads];
  // Calculate mean noise and mean signal
  cns_arg arg1[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg1[i].mask = mask;
    arg1[i].in1 = in1;
    arg1[i].in2 = in2;
    arg1[i].count = 0.0;
    arg1[i].noise = 0.0;
    arg1[i].power = 0.0;
    arg1[i].size = max;
    arg1[i].step = nthreads;
    arg1[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) calc_noise_signal_thread, &arg1[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  long double count = 0.0;
  long double noise = 0.0;
  long double power = 0.0;
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
    count += arg1[i].count;
    noise += arg1[i].noise;
    power += arg1[i].power;
  }
  noise /= count;
  power /= count;
  double psnr = fabsl(1.0 - noise / power);
  double rmsd = sqrtl(power);
  // Correct according to probability and power
  prob_arg arg2[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg2[i].mask = mask;
    arg2[i].in1 = in1;
    arg2[i].in2 = in2;
    arg2[i].out1 = out1;
    arg2[i].out2 = out2;
    arg2[i].res_stp_sd = node->stp / rmsd;
    arg2[i].size = max;
    arg2[i].step = nthreads;
    arg2[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) probability_correct_thread, &arg2[i])){
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
  node->max = psnr;
  return psnr;
}

void calc_noise_signal_thread(cns_arg *arg){
  int32_t i;
  double cur;
  for (i = arg->thread; i < arg->size; i += arg->step){
    // Normalise input transforms first
    arg->in1[i] = arg->in1[i] / arg->size;
    arg->in2[i] = arg->in2[i] / arg->size;
    // Do not calculate statistics from voxels outside the mask
    if (arg->mask->data[i] < 0.99){
      continue;
    }
    arg->count += 1.0;
    cur = arg->in1[i] - arg->in2[i];
    arg->noise += cur * cur;
    cur = arg->in1[i] + arg->in2[i];
    arg->power += cur * cur;
  }
  return;
}

void probability_correct_thread(prob_arg *arg){
  int32_t i;
  for (i = arg->thread; i < arg->size; i += arg->step){
    arg->out1[i] += arg->in1[i] * arg->res_stp_sd;
    arg->out2[i] += arg->in2[i] * arg->res_stp_sd;
  }
  return;
}
