
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

int get_num_jobs(void){
  // Obtain thread number from environmental variables
  char* thread_number = getenv("OMP_NUM_THREADS");
  int nthreads = 0;
  if (thread_number){
    // If thread number specified by user - use this one
    nthreads = atoi(thread_number);
  }
  if (nthreads < 1){
    // If thread number still not set - try sysconf
    nthreads = sysconf(_SC_NPROCESSORS_ONLN);
  }
  if (nthreads < 1){
    // If variables are both empty - use a single thread
    nthreads = 1;
  }
  return nthreads;
}

arguments *parse_args(int argc, char **argv){

  // Print usage and disclaimer

  printf("\n\n\n");

  char *splash =

    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n"
    "                                          #@#..*%@@@(                                                                                                          \n"
    "                                     .@#                .#@(                                    ,@@@(./#@@@(                                                   \n"
    "                                 %&                            &@*                          .@@               ,@%                                              \n"
    "                             *@,                                   @#                  ,@#  @                     @#                                           \n"
    "                          @%                                         .@.           .@/     @                         @                                         \n"
    "                         @                                             #(        @(       ,*                           #                                       \n"
    "                       .&                 *&*                     *@&    @    .@         /#                             %*                                     \n"
    "                      #,            (@,      &( /            (,.       /@*@   @         (/                                @                                    \n"
    "                      &    (    %@.(         *.,                      ,%   #*@         @.                                  *,                                  \n"
    "                      ( , @  %& (/          /%   *#&@@@@@&&&@@@@@%(    ( @  @         @                                     ,#                                 \n"
    "                     ,   % %/  .*,/                                 *&     @        %*                                       /,                                \n"
    "                     %   %@ /# .   ,/,,                           %       @        .                                          @                                \n"
    "                     @   @ .%#*@*                                @(     ,%                                                     @                               \n"
    "                     #  &   @.                                 @ &     %                                                        %                              \n"
    "                    /  ** @,      (,             .(         ,@   #    #                                                         ./                             \n"
    "                    @  @ @      /% *           @*&        %*   .,%@@@@(        #                                                 ,.                            \n"
    "                    / ,/&  @  ,,#%.   #%,..*@%% @       @   @( .//&,#@         @                                                  /                            \n"
    "                    @  %  %  @ .&  %/&%*..@* . @     /#  /.        @(%       & %                                                   **           &              \n"
    "                     @@   #*/  (,(      @*     @   .%               @/      .*.@&                                    @ &.            / .(&@@%/@#               \n"
    "                      &  %(           **      ./  &            , .( @       (   ##           @,    &                @  * @               ./@@,                 \n"
    "                     .,  @           @@&(#%@( *,@          @/@@&((@@@       @    @           #(    @ %,            %  #   @              /@                    \n"
    "                     %.  @         @& @     ,* @         %*    ,&   @       @     @         ,  @   (//& ,@@#.   (@&@  @   @                @&                  \n"
    "                     ,*           ,(* .(,  ##             @#  @*    @       @. #%/%@        *  (@*  &            , @   &  @                 @(                 \n"
    "                      &           &@       ,(#            #%@#,  .* (  /    &%%,%   @. /    # &( ,@  ,            ,%.  .# %                 ,@                 \n"
    "                    **@          @ @    ,        ,                  ,     (.@      %& *   #(  ,@,/&/            %%   . @                  .@                   \n"
    "                ,&    .&  *@    .*(&%/ % .      *&                   ( @ @ , ,&/   #   @*@  @   @# @/@@*           @    %  &.       #       (&                 \n"
    "                  ,@@@@  #@ . @    %   &      #*(  #(.             @/#  .@/@&@%   &   *&@ /@  @%  (*             &  %.@*   @      @#,     (,                   \n"
    "                    @   .@ % @  &   &   *   *                  .     ,@.  &  @/&&(/  %    #@%@# @&(   (,             @.*#@    @ %  .% (,   /@                  \n"
    "                   &*   @.(@/ /@@   @    ./ /      /% (%(,..  *@     .@/  @#   @*  &    ,    , .%@#(    %   #&       @,&.( &  / &  @   @   @                   \n"
    "                   @  %(   @,&( .,,  .(    &.      * ....,    #     @,#(% (.  #%.%% *@,            %      *,(@      (//(@ #    @/ @  , % .@,                   \n"
    "                   @ @       @@@   %   @    .@  &,         #&     &@#@   #@@  #* %       (         .      & @      ,&  /.*     @ @.   ,& //                    \n"
    "                   /&  .       (@*  #  /      */        .       @&&  %&      .    /@.                     @       *.  /&      %@* &      (&                    \n"
    "                  @,, #          &              (&  (%%*     (%#,&   % ,@(   ,     *   @(                 @      @   (&   /  .#  #      @,   ,@                \n"
    "             /(&(@     (      .%  ,.              @       (& .,.@     @.  @/  #   @     . @,              @    ,&   &    ,(      ,   ,@  *  ** ,.              \n"
    "           @    %   *        &  /.  %  *          @&   &@,     /       ,@   @ %  (,  @       @   ..       (,@@@%  .%     * ,    @   @        (@@%              \n"
    "          @     * (.  /     @,.  @   %            *,         @ @         @#   @#&(  /          @#          &     @      , .    ,  @*    (&       (@            \n"
    "         &      #  (   .   .,  %     **           ..        #.# %         *@   *@   %   *       ,@          ,&   @     / .     & @  /@*            (%          \n"
    "        #      %/%  (  @    %       @ (           *        ..  %@          //   @  (    (      %/#(       #@     %    #       @   #(          */     & .       \n"
    "        @,      , @ .%  ,   &   .,,/  @@         # @      .*  *.@           &   @. #    .    @,  .@     @/      %/           ,# &.        %.         #%        \n"
    "         @         @  &,@ ( %    @.(  (         @@  @    /.   / @            @  *,(        @*     @  &%         *@           & @      #@&*            @        \n"
    "         .(         % &(%**@    . (,   &       @@   @       /% /@            &    &      /@      ,@&            @,&          (@     .                 *,   .   \n"
    "          @             (@@      @@     .@      %   (#     %/ (*@            .&  /      .&    %( ,            @  @(         %@ @%,                     &       \n"
    "          %                      @@       &    &    (@    &,&  %@             @  @      @/    %*   & &      @    (.         , @,                       &       \n"
    "           ,         % .@*       @          @  @    @%@  &@    &@     ,@@@,  (@ .(          *  @  / %    %/       &         @                         @.       \n"
    "           @         .           #           *#&    @/*&@      @&   ,@,         @       .   (  @      ,@          @       #&                        /@         \n"
    "           @                    @              &.   &*( (      @   @(          @@           &  .                .  ,     @                         @/          \n"
    "           @(                  (*                ,@ /*( .     ,@, @           @ &           @ *,                   @   .@#@@&&,,                 &@ &          \n"
    "            .@@                @                     %* #     &%/@           & */           /@(.     &@(            % &                         @# @.          \n"
    "             .%               @                     %,  /     @%@           %. @       *     % (@                   %@*                   ,   %@* & (          \n"
    "              /.             @,                   /@          @#           @.  @       %     @   @                    @                   %  @@  &  @     .    \n"
    "               @           .@%                  /&@          ,(           @    @       %     .%  .%                   &                   # @*  @  #@     .    \n"
    "                *      ,@/  %                  % %@          &           @     @       #      @ @ @                   %                   &@   @  @.&          \n"
    "                @   @  .%  *                    & @    ,     @         /%      &       #        & (                   @                  @.   %  * &*          \n"
    "                 , ., #                        # #     &     @        &@@.     #       (       (,  %                  & @              .(.  &/ &   @           \n"
    "               / @ @  .  *                    /( @           @       @   @,    #       /       @   @                   %/            @#   %,( &   @*           \n"
    "          %@(    *  *&  .                     @ %,    .     ,#         *@,     /       *       @   @                               @(       %     @            \n"
    "          **      %     /                    %  @ /@%/.   %/        @@%     .,               &   /,                            @               %               \n"
    "           @(     @                         @  %.    #//    (/        / @     .*              /*    &                         .@                 @(            \n"
    "           %#      ,,                      @.  @            ,%          &%     &              &     @                       ##                   .@            \n"
    "            @        &&                  .&/  (/      * ,.   &           @,   .@              @     @                     @,                      (&           \n"
    "             %         .@                @%  ,@       %      @            @.  ,@ @@ @*        @     &*                #.@         *((/(//@,      @,(           \n"
    "             *#           @.            ,@   %&      .       &             @/ %@@   @*        @      @               %@.          ##,          .  .@           \n"
    "              (@            @,          @   ,/  #            (.               /    *@         @       *@            &%    (                 .     ,#           \n"
    "               %&             @*       ,%  &  .  #  *         &                    @@        #@          @         @                              %            \n"
    "                #@               *&@@@*   @  .     /          @                    @.       % @                %@&                              (%&            \n"
    "                 /@                     ,&//      #         *@@             (        .@@@&.   @                   .(                        *% . ,             \n"
    "                   @,                  @,    /@@&%/,*/#&@%*   *#          @(                  &,                   @                   /#                      \n"
    "                    @@               .@     &             .    @        #@                     @      .(%&*        &               .                           \n"
    "                      @@             / .@@, #            .&(&@ @      (@                        ((.   &    #        @                   .                      \n"
    "                        @@(        &         &/&/.  .*%%&*.         @@,                               %  .         &@                                          \n"
    "                          &@@@    &                               @(                             .    @@        .@@,                                           \n"
    "                             %@@@                              #&                                     ##  .,@&@@*                                              \n"
    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n"
    "                                            __   __   __    ____    __   ____   __    __  ______  ______  ____  _____                                          \n"
    "                                           /  /\\|  |\\|  \\  |   _|\\ /  /\\|    \\ |  |\\ |  ||      ||      ||   _||     \\                                         \n"
    "                                          /  / /|  |||   \\ |  |\\_\\/  / /|  >  ||  || |  ||__  __||_    _||  |\\ |  >  /\\                                        \n"
    "                                          \\  \\/ |  |||  , \\|   -| \\  \\/ |   _/\\|  || |  ||\\|  |\\_\\\\|  |\\\\|   -,|    / /                                        \n"
    "                                           \\  \\ |  |||  ` /|   -'\\ \\  \\ |  |\\\\/|  |_ |  || |  ||   |  || |   -\\|    \\/                                         \n"
    "                                           /  /\\|  |||   / |  |_\\/ /  /\\|  ||  |    ||  || |  ||   |  || |  |_\\|     \\                                         \n"
    "                                          /__/ /|__|||__/ /|____|\\/__/ /|__||  |____||__|| |__||   |__|| |____||__/\\__\\                                        \n"
    "                                          \\__\\/ \\__\\/\\__\\/ \\____\\/\\__\\/ \\__\\/  \\____\\\\__\\/ \\__\\/   \\__\\/ \\____\\/\\__\\\\__\\                                       \n"
    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n"
    "                                                                                                                                                               \n";


  printf("\n%s\n\n", splash);

  if (argc < 7){
    printf("\n    Usage: %s --v1 half_map1.mrc --v2 half_map2.mrc --o output_root --mask mask.mrc [ --spectrum || --rotfl ]\n\n", argv[0]);
  }

  printf("    PLEASE NOTE: SIDESPLITTER requires the unfiltered halfmaps and mask from each iteration or your results will be invalid\n");
  printf("                 Setting flag --spectrum outputs the natural SNR weighted spectrum rather than matching your input spectrum\n");
  printf("                 Setting flag --rotfl performs SNR tapering, matching input density in real-space rather than Fourier-space\n");
  printf("                 Remember - Junk in = Junk out! Please report any bug or observation to c.aylett@imperial.ac.uk, good luck!\n\n");
  printf("    SIDESPLITTER V1.2: LAFTER algorithm for halfmaps - 06-06-2020 GNU Public Licensed - K Ramlaul, CM Palmer and CHS Aylett\n\n");

  // Capture user requested settings
  int i;
  arguments *args = malloc(sizeof(arguments));
  memset(args, 0, sizeof(arguments));
  for (i = 1; i < argc; i++){
    if (!strcmp(argv[i], "--v1") && ((i + 1) < argc)){
      args->vol1 = argv[i + 1];
    } else if (!strcmp(argv[i], "--v2") && ((i + 1) < argc)){
      args->vol2 = argv[i + 1];
    } else if (!strcmp(argv[i], "--o") && ((i + 1) < argc)){
      args->out = argv[i + 1];
    } else if (!strcmp(argv[i], "--mask") && ((i + 1) < argc)){
      args->mask = argv[i + 1];
    } else if (!strcmp(argv[i], "--spectrum")){
      args->spec = 1;
    } else if (!strcmp(argv[i], "--rotfl")){
      args->rotf = 1;
    }
  }
  if (args->vol1 == NULL || args->vol2 == NULL){
    printf("    Necessary maps not found or unspecified - SIDESPLITTER absolutely requires the two halfset volumes and any mask applied\n\n");
    exit(1);
  }
  return args;
}

// Extend list by one using p-val
list *extend_list(list *current, double p){
  list *node = calloc(1, sizeof(list));
  node->res = current->res + current->stp;
  node->stp = p * (node->res / 64.0);
  node->prv = current;
  node->nxt = NULL;
  current->nxt = node;
  return node;
}

// End list for overfitting calculation
list *end_list(list *current){
  list *node = calloc(1, sizeof(list));
  node->res = current->res + current->stp;
  node->stp = 0.475 - (current->res + current->stp);
  if (node->stp < 0.0625){
    node->stp = 0.0625;
  }
  node->prv = current;
  node->nxt = NULL;
  current->nxt = node;
  return node;
}

// Read map header and data and return corresponding data structure
r_mrc *read_mrc(char *filename){

  FILE *f;
  int i;
  f = fopen(filename, "rb");
  if (!f){
    printf("\n\tError reading %s - bad file handle\n\n", filename);
    exit(1);
  }
  r_mrc *header = calloc(1, sizeof(r_mrc));
  fread(&header->n_crs, 4, 3, f);
  fread(&header->mode, 4, 1, f);
  /* We only accept mode 2 - c float32 / FORTRAN real - because implementing
     checks for other data types would require pointer casting everywhere */
  if (header->mode != 2){
    printf("Error reading %s - not 32 bit data \n", filename);
    exit(1);
  }
  fread(&header->start_crs,  4, 3,   f);
  fread(&header->n_xyz,      4, 3,   f);
  fread(&header->length_xyz, 4, 3,   f);
  fread(&header->angle_xyz,  4, 3,   f);
  fread(&header->map_crs,    4, 3,   f);
  fread(&header->d_min,      4, 1,   f);
  fread(&header->d_max,      4, 1,   f);
  fread(&header->d_mean,     4, 1,   f);
  fread(&header->ispg,       4, 1,   f);
  fread(&header->nsymbt,     4, 1,   f);
  fread(&header->extra,      4, 25,  f);
  fread(&header->ori_xyz,    4, 3,   f);
  fread(&header->map,        1, 4,   f);
  fread(&header->machst,     1, 4,   f);
  fread(&header->rms,        4, 1,   f);
  fread(&header->nlabl,      4, 1,   f);
  fread(&header->label,      1, 800, f);
  // Check cube
  if (header->n_crs[0] != header->n_crs[1] || header->n_crs[0] != header->n_crs[2]){
    printf("Error reading %s - map is not a cube\n", filename);
    exit(1);
  }
  /* Assign float array for data, initialize it to zero and read it in - note that the
     endianness is not corrected - architectures may therefore be cross-incompatible*/
  header->data = calloc((header->n_crs[0] * header->n_crs[1] * header->n_crs[2]), sizeof(float));
  if (!header->data){
    printf("Error reading %s - map not allocated\n", filename);
    exit(1);
  }
  fread(header->data, sizeof(float), (header->n_crs[0] * header->n_crs[1] * header->n_crs[2]), f);
  fclose(f);
  if (header->length_xyz[0] < 1e-9 || header->length_xyz[1] < 1e-9 || header->length_xyz[2] < 1e-9){
    header->length_xyz[0] = (float) header->n_xyz[0];
    header->length_xyz[1] = (float) header->n_xyz[1];
    header->length_xyz[2] = (float) header->n_xyz[2];
  }
  return header;
}

// Write MRC file given an mrc structure and corresponding data
void write_mrc(r_mrc* header, double *vol, char* filename, int32_t size){

  int i;
  double total    = size * size * size;
  double current  = 0;
  double tmp, sum = 0;

  float *hold = header->data;
  header->data = calloc(total, sizeof(float));
  if (!header->data){
    printf("Error reading %s - map not allocated\n", filename);
    exit(1);
  }

  for (i = 0; i < total; i++){
    // Convert double map to float for writing out
    header->data[i] = (float) vol[i];
  }

  // Calculate new min, max and mean figures for header
  header->d_min = header->data[0];
  header->d_max = header->data[0];
  for (i = 0; i < total; i++){
    current = (double) header->data[i];
    if (current < header->d_min){
      header->d_min = (float) current;
    }
    if (current > header->d_max){
      header->d_max = (float) current;
    }
    sum += current;
  }
  header->d_mean = (float) (sum / total);

  // Calculate RMSD to fill in the rms field for scaling
  sum = 0;
  for (i = 0; i < total; i++){
    tmp  = (double) (header->d_mean - header->data[i]);
    sum += tmp * tmp;
  }
  header->rms = (float) sqrt((sum / total));

  // Write out 1024 byte header
  FILE * f;
  f = fopen(filename, "wb");
  if (!f){
    printf("Error writing %s - bad file handle\n", filename);
    exit(1);
  }
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&header->mode,       4, 1,   f);
  fwrite(&header->start_crs,  4, 3,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&header->length_xyz, 4, 3,   f);
  fwrite(&header->angle_xyz,  4, 3,   f);
  fwrite(&header->map_crs,    4, 3,   f);
  fwrite(&header->d_min,      4, 1,   f);
  fwrite(&header->d_max,      4, 1,   f);
  fwrite(&header->d_mean,     4, 1,   f);
  fwrite(&header->ispg,       4, 1,   f);
  fwrite(&header->nsymbt,     4, 1,   f);
  fwrite(&header->extra,      4, 25,  f);
  fwrite(&header->ori_xyz,    4, 3,   f);
  fwrite(&header->map,        1, 4,   f);
  fwrite(&header->machst,     1, 4,   f);
  fwrite(&header->rms,        4, 1,   f);
  fwrite(&header->nlabl,      4, 1,   f);
  fwrite(&header->label,      1, 800, f);

  // Write data to file
  fwrite(header->data, 4, total, f);
  fclose(f);

  free(header->data);
  header->data = hold;
  return;
}
