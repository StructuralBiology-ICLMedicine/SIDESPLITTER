
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
  long double *spec1 = calloc(xyz, sizeof(long double));
  long double *spec2 = calloc(xyz, sizeof(long double));
  double maxres = get_spectrum(ki1, ki2, spec1, spec2, xyz, nthread);

  // Report FSC cut-off
  printf("\n\t FSC cut-off within mask = %12.6f \n", apix / maxres);

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Copy data into place
  add_map(vol1, ro1, nthread);
  add_map(vol2, ro2, nthread);

  // Execute forward transform
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Copy across ffts if tapering
  fftw_complex *inpk1 = NULL;
  fftw_complex *inpk2 = NULL;

  if (args->rotf){

    inpk1 = fftw_malloc(k_st);
    inpk2 = fftw_malloc(k_st);

    memcpy(inpk1, ki1, k_st);
    memcpy(inpk2, ki2, k_st);
  }
  
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

  double mean_p;

  // Noise suppression loop
  printf("\n\t Normalising -- Pass 1 \n");
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

    mean_p = normalise(ri1, ri2, ro1, ro2, mask, tail, xyz, nthread);
    
    if (tail->res + tail->stp >= maxres || mean_p <= 0.05){
      maxres = tail->res + tail->stp;
      break;
    }
          
    printf("\t Resolution = %12.6Lf | MeanProb = %12.6f | FSC = %12.6f \n", apix / (tail->res + tail->stp), mean_p, tail->fsc);
    fflush(stdout);

    tail = extend_list(tail, mean_p);

  } while (1);

  // Back-transform noise-suppressed maps
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Truncate by SNR
  printf("\n\t De-noising volume -- Pass 2 \n");
  printf("\n\t # Recovery indicates the fraction of the mask recovered by the current resolution");
  printf("\n\t # This should reach at least 1.0 but will preferably end up considerably higher -\n\n");
  fflush(stdout);

  char *name1 = "halfmap1.mrc";
  char *name2 = "halfmap2.mrc";
  
  if (args->out){
    size_t name_buffer = snprintf(NULL, 0, "%s%s", args->out, "_halfmap1.mrc") + 1;
    name1 = malloc(name_buffer);
    sprintf(name1, "%s%s", args->out, "_halfmap1.mrc");
    name2 = malloc(name_buffer);
    sprintf(name2, "%s%s", args->out, "_halfmap2.mrc");
  }

  // Choose tapering loop if required
  if (args->rotf){

    double *ori1 = fftw_malloc(r_st);
    double *ori2 = fftw_malloc(r_st);

    memset(ori1, 0, r_st);
    memset(ori2, 0, r_st);

    fftw_complex *oki1 = fftw_malloc(k_st);
    fftw_complex *oki2 = fftw_malloc(k_st);

    memset(oki1, 0, k_st);
    memset(oki2, 0, k_st);

    fftw_plan fft_oki1_ori1 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, oki1, ori1, FFTW_ESTIMATE);
    fftw_plan fft_oki2_ori2 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, oki2, ori2, FFTW_ESTIMATE);
  
    do {

      lowpass_filter(ki1, ko1, tail, xyz, nthread);
      lowpass_filter(ki2, ko2, tail, xyz, nthread);

      lowpass_filter(inpk1, oki1, tail, xyz, nthread);
      lowpass_filter(inpk2, oki2, tail, xyz, nthread);

      fftw_execute(fft_ko1_ri1);
      fftw_execute(fft_ko2_ri2);

      fftw_execute(fft_oki1_ori1);
      fftw_execute(fft_oki2_ori2);

      mean_p = taper_map(ri1, ri2, ro1, ro2, ori1, ori2, mask, tail, args, xyz, nthread);

      printf("\t Resolution = %12.6Lf | Recovery = %12.6f\n", apix / (tail->res + tail->stp), mean_p);
      fflush(stdout);
    
      if (tail->prv == NULL){
	break;
      } else{
	tail = tail->prv;
      }

    } while (1);

    // Renormalise maps
    int32_t total = xyz * xyz * xyz;

    for (i = 0; i < total; i++){
      ro1[i] /= (double) total;
      ro2[i] /= (double) total;
    }

    // Output maps if SNR tapering
    write_mrc(vol1, ro1, name1, xyz);
    write_mrc(vol2, ro2, name2, xyz);

    // Over and out...
    printf("\n\n\n\t ++++ ++++ That's All Folks! ++++ ++++ \n\n\n");

    return 0;
    
  } else {
    do {

      lowpass_filter(ki1, ko1, tail, xyz, nthread);
      lowpass_filter(ki2, ko2, tail, xyz, nthread);

      fftw_execute(fft_ko1_ri1);
      fftw_execute(fft_ko2_ri2);

      mean_p = truncate_map(ri1, ri2, ro1, ro2, mask, tail, args, xyz, nthread);

      printf("\t Resolution = %12.6Lf | Recovery = %12.6f\n", apix / (tail->res + tail->stp), mean_p);
      fflush(stdout);
    
      if (tail->prv == NULL){
	break;
      } else{
	tail = tail->prv;
      }

    } while (1);
  }

  // Back-transform noise-suppressed maps
  fftw_execute(fft_ro1_ki1);
  fftw_execute(fft_ro2_ki2);

  // Zero fill maps
  memset(ro1, 0, r_st);
  memset(ro2, 0, r_st);

  // Noise suppression loop 2
  printf("\n\t Reapplying spectum \n");
  printf("\n\t # Spectrum indicates the spectral power reapplied at the current resolution\n\n");
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

    fftw_execute(fft_ko1_ri1);
    fftw_execute(fft_ko2_ri2);

    reverse_norm(ri1, ri2, ro1, ro2, mask, tail, xyz, nthread);

    printf("\t Resolution = %12.6Lf | Spectrum = %12.6Lf \n", apix / (tail->res + tail->stp), tail->pwr);
    fflush(stdout);

    if (tail->nxt == NULL){
      break;
    } else{
      tail = tail->nxt;
    }
    
  } while (1);

  // Apply masks in situ
  apply_mask(mask, ro1, nthread);
  apply_mask(mask, ro2, nthread);

  // Output final volume
  printf("\n\t Writing noise truncated MRC files\n");
  fflush(stdout);
  
  if (!args->spec){

    fftw_execute(fft_ro1_ki1);
    fftw_execute(fft_ro2_ki2);

    apply_spectrum(ki1, ki2, spec1, spec2, maxres, xyz, nthread);

    fftw_plan fft_ki1_ri1 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ki1, ri1, FFTW_ESTIMATE);
    fftw_plan fft_ki2_ri2 = fftw_plan_dft_c2r_3d(xyz, xyz, xyz, ki2, ri2, FFTW_ESTIMATE);

    fftw_execute(fft_ki1_ri1);
    fftw_execute(fft_ki2_ri2);
    
    // Renormalise maps
    int32_t total = xyz * xyz * xyz;

    for (i = 0; i < total; i++){
      ri1[i] /= (double) total;
      ri2[i] /= (double) total;
    }

    write_mrc(vol1, ri1, name1, xyz);
    write_mrc(vol2, ri2, name2, xyz);

  } else {

    // Renormalise maps
    int32_t total = xyz * xyz * xyz;

    for (i = 0; i < total; i++){
      ro1[i] /= (double) total;
      ro2[i] /= (double) total;
    }

    write_mrc(vol1, ro1, name1, xyz);
    write_mrc(vol2, ro2, name2, xyz);

  }

  // Over and out...
  printf("\n\n\n\t ++++ ++++ That's All Folks! ++++ ++++ \n\n\n");

  return 0;
}
