#!/bin/bash

# Wrapper to run relion_external_reconstruct and sidesplitter
# Version 1.0
# Author: Colin M. Palmer

# Usage:
#     sidesplitter_wrapper.sh path/to/relion_external_reconstruct_star_file.star
#
# To use from RELION 3.1, set the RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE environment variable to point to this script,
# set SIDESPLITTER to point to the sidesplitter binary (or make sure sidesplitter can be found via PATH), and run
# relion_refine with the --external_reconstruct argument. For example:
#
#     export RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE=/path/to/sidesplitter_wrapper.sh
#     export SIDESPLITTER=/path/to/sidesplitter
#
# then run RELION auto-refine from the GUI and put "--external_reconstruct" in the additional arguments box. To run on
# a cluster, depending on your configuration you might need to put the environment variable definitions into your
# submission script.

# Troubleshooting
#
# If you have problems running SIDESPLITTER using this script, the first thing to check is that external reconstruction
# from RELION is working correctly. Try running a normal refinement job, using the "--external_reconstruct" argument
# but without setting the RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE environment variable. If this fails, the problem is
# likely to be with your RELION installation - perhaps it is the wrong version, or different installations are
# conflicting with each other. If normal external reconstruction is successful, the problem is likely to be with the
# SIDESPLITTER installation, or a bug in this script.

# How this script works:
#
# If the target file name contains "_half", this script assumes two copies of itself will be running (for the two half
# data sets). The two copies will coordinate with each other by creating, checking for and deleting a directory called
# "sidesplitter_running". Both scripts will run relion_external_reconstruct for their given half data set. The first
# script will wait for both reconstructions to finish, then call SIDESPLITTER to process both half maps. The second
# script will wait for the first to finish running SIDESPLITTER and then exit (because if either of the scripts exits
# before the processing is finished, RELION moves on and tries to continue its own processing before the filtered
# volumes are ready).
#
# If the target file name does not contain "_half", this script assumes there is only a single copy of itself running.
# In this case it call relion_external_reconstruct, waits for the reconstruction to finish and then exits.
# This handles the final iteration when the two half sets are combined, at which point RELION calls the external
# reconstruction program just once to reconstruct the final combined volume.
#
# Note that this script is not particularly robust. If one of the commands fails, it's possible you might need to
# manually tidy up the job directory and remove the "sidesplitter_running" directory to avoid problems in the next run.


#### Configuration

# SIDESPLITTER command ("sidesplitter" unless defined in the SIDESPLITTER environment variable)
sidesplitter_default=sidesplitter
sidesplitter=${SIDESPLITTER:-${sidesplitter_default}}

# Change to true to activate debug output (to check for problems with process coordination)
debug=false


#### Main script

# Expect this to be called with one argument, pointing at a STAR file for relion_external_reconstruct
base_path=${1%.star}

running_ind_dir=${base_path%/*}/sidesplitter_running

$debug && echo "$$ Checking for existence of $running_ind_dir ..."

if mkdir "$running_ind_dir" 2> /dev/null; then
  first=true
  $debug && echo "$$ Created $running_ind_dir"
else
  first=false
  $debug && echo "$$ $running_ind_dir already exists"
fi

echo "$$ $(date) Running relion_external_reconstruct $1 > $base_path.out 2> $base_path.err"
relion_external_reconstruct "$1" > "$base_path.out" 2> "$base_path.err"

if [[ $base_path != *"_half"* ]]; then
  if $first; then
    $debug && echo "$$ $(date) Only a single reconstruction, removing $running_ind_dir and exiting."
    rmdir "$running_ind_dir"
    exit 0
  else
    $debug && echo "$$ $(date) Error! Found pre-existing $running_ind_dir for single reconstruction job."
    exit 1
  fi
fi

$debug && echo "$$ Moving output file ${base_path}.mrc to ${base_path}_orig.mrc"
mv "${base_path}.mrc" "${base_path}_orig.mrc"

if $first; then
  $debug && echo "$$ $(date) First reconstruct job finished; waiting for $running_ind_dir to disappear"
  while [[ -d $running_ind_dir ]]; do
    $debug && echo "$$ $(date) $running_ind_dir still exists; waiting..."
    sleep 5
  done

  $debug && echo "$$ $(date) $running_ind_dir has disappeared. Moving on to sidesplitter step."

  mask=$(awk '/fn_mask/ { gsub(/^"|"$/, "", $2) ; print $2 }' ${base_path%/*}/job.star)

  if [[ -z "$mask" ]]; then
    echo "Warning: no mask found! SIDESPLITTER will give better results if you use a mask."
    echo "$$ $(date) Running $sidesplitter  --v1 ${base_path/half2/half1}_orig.mrc --v2 ${base_path/half1/half2}_orig.mrc > ${base_path%_half*}_sidesplitter.out"
    $sidesplitter  --v1 "${base_path/half2/half1}_orig.mrc" --v2 "${base_path/half1/half2}_orig.mrc" > "${base_path%_half*}_sidesplitter.out"
  else
    echo "$$ $(date) Running $sidesplitter  --v1 ${base_path/half2/half1}_orig.mrc --v2 ${base_path/half1/half2}_orig.mrc --mask $mask > ${base_path%_half*}_sidesplitter.out"
    $sidesplitter  --v1 "${base_path/half2/half1}_orig.mrc" --v2 "${base_path/half1/half2}_orig.mrc" --mask "$mask" > "${base_path%_half*}_sidesplitter.out"
  fi

  $debug && echo "$$ Moving sidesplitter output halfmap1.mrc to ${base_path/half2/half1}.mrc"
  mv halfmap1.mrc "${base_path/half2/half1}.mrc"

  $debug && echo "$$ Moving sidesplitter output halfmap2.mrc to ${base_path/half1/half2}.mrc"
  mv halfmap2.mrc "${base_path/half1/half2}.mrc"

  $debug && echo "$$ $(date) Finished sidesplitter. Recreating $running_ind_dir to signal job finished."
  mkdir "$running_ind_dir"

else
  $debug && echo "$$ $(date) Second reconstruct job finished; removing $running_ind_dir"
  rmdir "$running_ind_dir"
  $debug && echo "$$ $(date) Waiting for $running_ind_dir to reappear"
  while [[ ! -d $running_ind_dir ]]; do
    $debug && echo "$$ $(date) $running_ind_dir does not exist; waiting..."
    sleep 60
  done
  $debug && echo "$$ $(date) $running_ind_dir has reappeared. Removing it and exiting"
  rmdir "$running_ind_dir"
fi

