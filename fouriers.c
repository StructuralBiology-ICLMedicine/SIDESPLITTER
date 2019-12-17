
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
#include "fouriers.h"

// Add FFT in to FFT out
void add_fft(fftw_complex *in, fftw_complex *out, int32_t full, int32_t nthreads){
  int32_t size = full * full * (full / 2 + 1), i;
  pthread_t threads[nthreads];
  add_fft_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) add_fft_thread, &arg[i])){
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

void add_fft_thread(add_fft_arg *arg){
  for(int32_t index = arg->thread; index < arg->size; index += arg->step){
    arg->out[index] += arg->in[index];
  }
  return;
}

// Calculate FSC over map
double calc_fsc(fftw_complex *half1, fftw_complex *half2, int32_t full, int32_t nthreads){
  int32_t size = full * full * (full / 2 + 1), i;
  pthread_t threads[nthreads];
  calc_fsc_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].half1 = half1;
    arg[i].half2 = half2;
    arg[i].numerator = 0.0;
    arg[i].denomin_2 = 0.0;
    arg[i].denomin_1 = 0.0;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) calc_fsc_thread, &arg[i])){
      printf("\nThread initialisation failed!\n");
      fflush(stdout);
      exit(1);
    }
  }
  long double numerator = 0.0;
  long double denomin_1 = 0.0;
  long double denomin_2 = 0.0;
  // Join threads
  for (i = 0; i < nthreads; i++){
    if (pthread_join(threads[i], NULL)){
      printf("\nThread failed during run!\n");
      fflush(stdout);
      exit(1);
    }
    numerator += arg[i].numerator;
    denomin_1 += arg[i].denomin_1;
    denomin_2 += arg[i].denomin_2;
  }
  return (double) (numerator / sqrtl(fabsl(denomin_1 * denomin_2)));
}

void calc_fsc_thread(calc_fsc_arg *arg){
  for(int32_t index = arg->thread; index < arg->size; index += arg->step){
    arg->numerator += creal(arg->half1[index] * conj(arg->half2[index]));
    arg->denomin_1 += creal(arg->half1[index] * conj(arg->half1[index]));
    arg->denomin_2 += creal(arg->half2[index] * conj(arg->half2[index]));
  }
  return;
}

// Apply bandpass to in and writes to out
void bandpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t full, int32_t nthreads){
  double hires = node->res + node->stp;
  double lores = node->res;
  hires = hires * hires;
  lores = lores * lores;
  double dim = (double) full;
  int32_t size = (full / 2) + 1, i;
  int32_t full_size = full * size;
  pthread_t threads[nthreads];
  filter_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].hires = hires;
    arg[i].lores = lores;
    arg[i].dim = dim;
    arg[i].full_size = full_size;
    arg[i].full = full;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) bandpass_filter_thread, &arg[i])){
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

void bandpass_filter_thread(filter_arg* arg){
  double norms, kd, jd, id;
  int32_t index;
  for(int _k = 0, k = 0; _k < arg->full; _k++, k = (_k < arg->size) ? _k : _k - arg->full){
    kd = ((double) k) / arg->dim;
    for(int _j = 0, j = 0; _j < arg->full; _j++, j = (_j < arg->size) ? _j : _j - arg->full){
      jd = ((double) j) / arg->dim;
      for(int _i = arg->thread, i = arg->thread; _i < arg->size; _i += arg->step, i = _i){
        id = ((double) i) / arg->dim;
        norms = kd * kd + jd * jd + id * id;
        index = _k * arg->full_size + _j * arg->size + _i;
        arg->out[index] = arg->in[index] * (sqrt(1.0 / (1.0 + pow((norms / arg->hires), 8.0))) - sqrt(1.0 / (1.0 + pow((norms / arg->lores), 8.0))));
      }
    }
  }
  return;
}

// Butterworth lowpass from in to out
void lowpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t full, int32_t nthreads){
  double hires = node->res + node->stp;
  hires = hires * hires;
  double dim = (double) full;
  int32_t size = (full / 2) + 1, i;
  int32_t full_size = full * size;
  pthread_t threads[nthreads];
  filter_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in = in;
    arg[i].out = out;
    arg[i].hires = hires;
    arg[i].dim = dim;
    arg[i].full_size = full_size;
    arg[i].full = full;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) lowpass_filter_thread, &arg[i])){
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

void lowpass_filter_thread(filter_arg *arg){
  double norms, kd, jd, id;
  int32_t index;
  for(int _k = 0, k = 0; _k < arg->full; _k++, k = (_k < arg->size) ? _k : _k - arg->full){
    kd = ((double) k) / arg->dim;
    for(int _j = 0, j = 0; _j < arg->full; _j++, j = (_j < arg->size) ? _j : _j - arg->full){
      jd = ((double) j) / arg->dim;
      for(int _i = arg->thread, i = arg->thread; _i < arg->size; _i += arg->step, i = _i){
        id = ((double) i) / arg->dim;
        norms = kd * kd + jd * jd + id * id;
        index = _k * arg->full_size + _j * arg->size + _i;
        arg->out[index] = arg->in[index] * sqrt(1.0 / (1.0 + pow((norms / arg->hires), 8.0)));
      }
    }
  }
  return;
}

// Calculate spectrum over map
double get_spectrum(fftw_complex *half1, fftw_complex *half2, double *spec1, double *spec2, int32_t full, int32_t nthreads){
  double cut = 0.0;
  double dim = (double) full;
  int32_t size = (full / 2) + 1, i, j;
  int32_t full_size = full * size;
  int32_t *n = calloc(full, sizeof(int32_t));
  double *sum = calloc(full, sizeof(double));
  double *sub = calloc(full, sizeof(double));
  pthread_t threads[nthreads];
  spec_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in1 = half1;
    arg[i].in2 = half2;
    arg[i].out1 = calloc(full, sizeof(double));
    arg[i].out2 = calloc(full, sizeof(double));
    arg[i].n = calloc(full, sizeof(int32_t));
    arg[i].sum = calloc(full, sizeof(double));
    arg[i].sub = calloc(full, sizeof(double));
    arg[i].dim = dim;
    arg[i].full_size = full_size;
    arg[i].full = full;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) get_spec_thread, &arg[i])){
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
    for (j = 0; j < full; j++){
      n[j] += arg[i].n[j];
      spec1[j] += arg[i].out1[j];
      spec2[j] += arg[i].out2[j];
      sum[j] += arg[i].sum[j];
      sub[j] += arg[i].sub[j];
    }
    free(arg[i].n);
    free(arg[i].out1);
    free(arg[i].out2);
    free(arg[i].sum);
    free(arg[i].sub);
  }
  for (i = 0; i < full; i++){
    spec1[i] = (spec1[i] / (double) n[i]);
    spec2[i] = (spec2[i] / (double) n[i]);
    if ((cut > 0.0) || (spec1[i] > 0.0 && spec2[i] > 0.0 && spec1[i] < 0.1 && spec2[i] < 0.1) || (log2(sum[i] / sub[i]) < 0.25)){
      if (cut == 0.0){
	cut = ((double) i) / (full * 2.0);
      }
      spec1[i] = 0.0;
      spec2[i] = 0.0;
      continue;
    }
  }
  free(n);
  if (cut > 0.0){
    return cut;
  }
  return 0.45;
}

void get_spec_thread(spec_arg *arg){
  double kd, jd, id;
  int32_t index, norms;
  for(int _k = 0, k = 0; _k < arg->full; _k++, k = (_k < arg->size) ? _k : _k - arg->full){
    kd = (double) k;
    for(int _j = 0, j = 0; _j < arg->full; _j++, j = (_j < arg->size) ? _j : _j - arg->full){
      jd = (double) j;
      for(int _i = arg->thread, i = arg->thread; _i < arg->size; _i += arg->step, i = _i){
        id = (double) i;
        norms = (int32_t) (sqrt(fabs(kd * kd + jd * jd + id * id)) * 2.0);
	if (norms >= arg->full){
	  continue;
	}
        index = _k * arg->full_size + _j * arg->size + _i;
        arg->out1[norms] += sqrt(fabs(creal(arg->in1[index] * conj(arg->in1[index]))));
	arg->out2[norms] += sqrt(fabs(creal(arg->in2[index] * conj(arg->in2[index]))));
	if (arg->sum && arg->sub){
	  arg->sum[norms] += creal((arg->in1[index] + arg->in2[index]) * conj(arg->in1[index] + arg->in2[index]));
	  arg->sub[norms] += creal((arg->in1[index] - arg->in2[index]) * conj(arg->in1[index] - arg->in2[index]));
	}
	arg->n[norms]++;
      }
    }
  }
  return;
}

// Apply spectrum over map
void apply_spectrum(fftw_complex *half1, fftw_complex *half2, double *spec1, double *spec2, double maxres, int32_t full, int32_t nthreads){
  double dim = (double) full;
  int32_t size = (full / 2) + 1, i, j;
  int32_t full_size = full * size;
  int32_t *n = calloc(full, sizeof(int32_t));
  double *cor1 = calloc(full, sizeof(double));
  double *cor2 = calloc(full, sizeof(double));
  pthread_t threads[nthreads];
  spec_arg arg[nthreads];
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].in1 = half1;
    arg[i].in2 = half2;
    arg[i].out1 = calloc(full, sizeof(double));
    arg[i].out2 = calloc(full, sizeof(double));
    arg[i].n = calloc(full, sizeof(int32_t));
    arg[i].sum = NULL;
    arg[i].sub = NULL;
    arg[i].dim = dim;
    arg[i].full_size = full_size;
    arg[i].full = full;
    arg[i].size = size;
    arg[i].step = nthreads;
    arg[i].thread = i;
    if (pthread_create(&threads[i], NULL, (void*) get_spec_thread, &arg[i])){
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
    for (j = 0; j < full; j++){
      n[j] += arg[i].n[j];
      cor1[j] += arg[i].out1[j];
      cor2[j] += arg[i].out2[j];
    }
    free(arg[i].n);
    free(arg[i].out1);
    free(arg[i].out2);
  }
  // Normalise
  int32_t cut = (int32_t) (maxres * full * 2.0);
  for (i = 0; i < full; i++){
    if (i < cut){
      cor1[i] = spec1[i] / (cor1[i] / (double) n[i]);
      cor2[i] = spec2[i] / (cor2[i] / (double) n[i]);
    } else {
      cor1[i] = 0.0;
      cor2[i] = 0.0;
    }
  }
  // Start threads
  for (i = 0; i < nthreads; i++){
    arg[i].out1 = cor1;
    arg[i].out2 = cor2;
    if (pthread_create(&threads[i], NULL, (void*) apply_spec_thread, &arg[i])){
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
  free(cor1);
  free(cor2);
  free(n);
  return;
}

void apply_spec_thread(spec_arg *arg){
  double kd, jd, id;
  int32_t index, norms;
  for(int _k = 0, k = 0; _k < arg->full; _k++, k = (_k < arg->size) ? _k : _k - arg->full){
    kd = (double) k;
    for(int _j = 0, j = 0; _j < arg->full; _j++, j = (_j < arg->size) ? _j : _j - arg->full){
      jd = (double) j;
      for(int _i = arg->thread, i = arg->thread; _i < arg->size; _i += arg->step, i = _i){
        id = (double) i;
        norms = (int32_t) (sqrt(fabs(kd * kd + jd * jd + id * id)) * 2.0);
        index = _k * arg->full_size + _j * arg->size + _i;
	if (norms >= arg->full){
	  arg->in1[index] *= 0.0;
	  arg->in2[index] *= 0.0;
	  continue;
	}
        arg->in1[index] *= arg->out1[norms];
	arg->in2[index] *= arg->out2[norms];
      }
    }
  }
  return;
}
