
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
 * Program: SIDESPLITTER V1.0                                               
 *                                                                         
 * Authors: Chris Aylett                                                   
 *          Colin Palmer                                                   
 *                                                                         
 */

// Library header inclusion for linking                                  
#include "sidesplitter.h"

// Main algorithm function
int main(int argc, char **argv){

  int32_t i, j;

  // Get arguments
  arguments *args = parse_args(argc, argv);
  int32_t nthread = get_num_jobs();
  
  // Read MRC inputs
  r_mrc *vol1 = read_mrc(args->vol1);
  r_mrc *vol2 = read_mrc(args->vol2);
  r_mrc *mask;
  if (args->mask){
    mask = read_mrc(args->mask);
  } else {
    mask = make_msk(vol1, (double) vol1->n_crs[0] / 4, nthread);
  }
  double apix = vol1->length_xyz[0] / (float) vol1->n_xyz[0];

  // Check map sizes and CPUs
  int32_t xyz = mask->n_crs[0];
  if (xyz != vol1->n_crs[0] || xyz != vol2->n_crs[0]){
    printf("\n\t MAPS MUST BE THE SAME SIZE! \n");
    return 1;
  }

  size_t r_st = xyz * xyz * xyz * sizeof(double);
  size_t k_st = xyz * xyz * ((xyz / 2) + 1) * sizeof(fftw_complex);

  // FFTW set-up
  printf("\n\t Setting up threads and maps\n");
  printf("\n\t Using %i threads. If you want to override this, set the OMP_NUM_THREADS environment variable.\n", nthread);
  fftw_init_threads();
  fftw_plan_with_nthreads(nthread);

  // Allocate memory for maps
  double *ri1 = fftw_malloc(r_st);
  double *ri2 = fftw_malloc(r_st);
  double *ro1 = fftw_malloc(r_st);
  double *ro2 = fftw_malloc(r_st);
  fftw_complex *ki1 = fftw_malloc(k_st);
  fftw_complex *ki2 = fftw_malloc(k_st);
  fftw_complex *ko1 = fftw_malloc(k_st);
  fftw_complex *ko2 = fftw_malloc(k_st);
  
  // Make FFTW plans
  printf("\n\t FFTW doing its thing - ");
  fflush(stdout);
  fftw_plan fft_ro1_ki1 = fftw_plan_dft_r2c_3d(xyz, xyz, xyz, ro1, ki1, FFTW_MEASURE);
  printf("#");
  fflush(stdout);
  fftw_plan fft_ro2_ki2 = fftw_plan_dft_r2c_3d(xyz, xyz, xyz, ro2, ki2, FFTW_ESTIMATE);
  printf("#");
  fflush(stdout);
  fftw_plan fft_ko1_ri1 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ko1, ri1, FFTW_MEASURE);
  printf("#");
  fflush(stdout);
  fftw_plan fft_ko2_ri2 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ko2, ri2, FFTW_ESTIMATE);
  printf("#\n");
  fflush(stdout);

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);
  memset(ri1, 0, r_st);
  memset(ri2, 0, r_st);
  memset(ko1, 0, k_st);
  memset(ko2, 0, k_st);
  memset(ki1, 0, k_st);
  memset(ki2, 0, k_st);

  // Copy data into place
  add_map(vol1, ro1, nthread);
  add_map(vol2, ro2, nthread);

  // Apply masks in situ
  apply_mask(mask, ro1, nthread);
  apply_mask(mask, ro2, nthread);
  
  // Execute forward transform
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Obtain spectra
  double *spec1 = calloc(xyz, sizeof(double));
  double *spec2 = calloc(xyz, sizeof(double));
  double maxres = get_spectrum(ki1, ki2, spec1, spec2, xyz, nthread);

  // Report FSC cut-off
  printf("\n\t FSC = 0.143 within mask = %12.6f \n", apix / maxres);

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Copy data into place
  add_map(vol1, ro1, nthread);
  add_map(vol2, ro2, nthread);

  // Execute forward transform
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Zero centre
  ki1[0] = 0.0 + 0.0I;
  ki2[0] = 0.0 + 0.0I;

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Initialise list
  list head;
  head.res = 0.000;
  head.stp = 0.025;
  head.prv = NULL;
  head.nxt = NULL;
  head.crf = 0.00;
  head.fsc = 1.00;

  list *tail = &head;

  // Loop parameterisation
  double prev_p = 1.00;
  double mean_p = 1.00;
  double step_p = 0.00;

  // Noise suppression loop
  printf("\n\t Suppressing noise -- Pass 1 \n");
  printf("\n\t # Resolution is reported in Ångströms [Å] everywhere it is quoted ");
  printf("\n\t # MeanProb records the estimated probability voxels are not noise ");
  printf("\n\t # FSC indicates the Fourier Shell Correlation between half sets -\n\n");
  fflush(stdout);

  i = 0;
  do {
    if (tail->res == 0.0){
      lowpass_filter(ki1, ko1, tail, xyz, nthread);
      lowpass_filter(ki2, ko2, tail, xyz, nthread);
    } else {
      bandpass_filter(ki1, ko1, tail, xyz, nthread);
      bandpass_filter(ki2, ko2, tail, xyz, nthread);
    }

    tail->fsc = calc_fsc(ko1, ko2, xyz, nthread);
    tail->crf = sqrt(fabs((2.0 * tail->fsc) / (1.0 + tail->fsc)));

    fftw_execute(fft_ko1_ri1);
    fftw_execute(fft_ko2_ri2);

    mean_p = suppress_noise(ri1, ri2, ro1, ro2, mask, tail, xyz, nthread);
    
    if (tail->res + tail->stp > 0.475 || (mean_p < 0.0 && tail->res + tail->stp > maxres)){
      maxres = tail->res + tail->stp;
      break;
    } else if (mean_p < 0.0){
      mean_p *= -1.0;
    }
    
    printf("\t Resolution = %12.6f | MeanProb = %12.6f | FSC = %12.6f \n", apix / (tail->res + tail->stp), mean_p, tail->fsc);
    fflush(stdout);

    tail = extend_list(tail, mean_p);
    
    prev_p = mean_p;

  } while (1);

  // Back-transform noise-suppressed maps
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Zero centre
  ki1[0] = 0.0 + 0.0I;
  ki2[0] = 0.0 + 0.0I;

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Truncate by SNR
  printf("\n\t De-noising volume -- Pass 2 \n");
  printf("\n\t # Recovery indicates the fraction of the mask recovered by the current resolution");
  printf("\n\t # This should reach at least 1.0 but will preferably end up considerably higher -\n\n");
  fflush(stdout);
  do {

    lowpass_filter(ki1, ko1, tail, xyz, nthread);
    lowpass_filter(ki2, ko2, tail, xyz, nthread);

    fftw_execute(fft_ko1_ri1);
    fftw_execute(fft_ko2_ri2);

    mean_p = truncate_map(ri1, ri2, ro1, ro2, mask, tail, args, xyz, nthread);

    printf("\t Resolution = %12.6f | Recovery = %12.6f\n", apix / (tail->res + tail->stp), mean_p);
    fflush(stdout);

    if (tail->prv == NULL){
      break;
    } else{
      tail = tail->prv;
    }

  } while (1);

  // Apply masks in situ
  apply_mask(mask, ro1, nthread);
  apply_mask(mask, ro2, nthread);
  
  // Forward transform
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  apply_spectrum(ki1, ki2, spec1, spec2, maxres, xyz, nthread);

  // Output final volume
  printf("\n\t Writing noise truncated MRC files\n");
  fflush(stdout);

  char *name1 = "halfmap1.mrc";
  char *name2 = "halfmap2.mrc";

  fftw_plan fft_ki1_ri1 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ki1, ri1, FFTW_ESTIMATE);
  fftw_plan fft_ki2_ri2 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ki2, ri2, FFTW_ESTIMATE);

  fftw_execute(fft_ki1_ri1);
  fftw_execute(fft_ki2_ri2);

  write_mrc(vol1, ri1, name1, xyz);
  write_mrc(vol2, ri2, name2, xyz);

  // Over and out...
  printf("\n\n\n\t ++++ ++++ That's All Folks! ++++ ++++ \n\n\n");

  return 0;
}
