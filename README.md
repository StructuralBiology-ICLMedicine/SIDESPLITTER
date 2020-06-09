
# SIDESPLITTER [ Local SNR Filter for both SIDES of SPLIT refinements ]
KAILASH RAMLAUL, COLIN PALMER, TAKANORI NAKANE & CHRISTOPHER H S AYLETT
09/06/2020


## Important points:
- SIDESPLITTER is a local SNR filter for Single Particle reconstruction
- SIDESPLITTER is intended to maintain independence between refinements
- SIDESPLITTER requires independent half-maps and any binary mask used
- SIDESPLITTER outputs denoised half volumes for further alignement
- SIDESPLITTER requires FFTW3 with POSIX threads for compilation
- SIDESPLITTER is a C99 program (a compilation script is included)
- SIDESPLITTER prints all its available options to stdout when called


## Usage of SIDESPLITTER with RELION
- A wrapper script is provided with SIDESPLITTER for ease of RELION use
- The estimation of the improvment in SNR is performed by the wrapper
- This improvement is estimated according to modification of the FSC:

  FSC_SS = (sqrt(2 (FSC + FSC^2)) + FSC) / (2 + FSC)

- Because SNR improvement is local this can be an over or underestimate
- If RELION refinement becomes unstable it may help to remove this


## Testing and Feedback
- SIDESPLITTER is a relatively novel technique. We believe that further
  validation of our method will best be facilitated by widespread use,
  and would actively encourage users to communicate any results from
  particularly difficult or interesting SIDESPLITTER refinements,
  especially heterogenous, conformationally flexible, poorly alignable,
  or contaminated datasets.
- Novel issues, interesting outcomes, or bugs in SIDESPLITTER should be
  reported to the authors or CCPEM. E-mail c.aylett@imperial.ac.uk .


## Installation and License
- SIDESPLITTER is a performance optimised C program using FFTW3 Fourier
  transformation (Frigo and Johnson, 2005) for speed and portability
- Compilation is easiest with CMake. Build and run SIDESPLITTER with the
  commands:

```bash
mkdir build && cd build
cmake ..
make
./sidesplitter
```

- SIDESPLITTER is open source and is made available under the GNU public
  license, which should be included in any package.

  ++ Good luck & happy cryo-EM-ing - Kailash, Colin, Takanori & Chris ++
