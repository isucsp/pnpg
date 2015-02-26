## Image Reconstruction Source Code

### What does this package contains

1. Accelerated Proximal-Gradient Algorithms

Reconstruction of Nonnegative Sparse Signals Using Accelerated Proximal-Gradient Algorithms

Nesterov's Proximal-Gradient Algorithms for Reconstructing Nonnegative Signals with Sparse Transform Coefficients

2. Beam Hardening Correction Algorithms


### How to Install

To install this package, first download the repository by running

    git clone https://github.com/isucsp/imgRecSrc.git

after downloading, from MATLAB change your current folder to `imgRecSrc/` and
execute `setupPath.m` file to add the paths to the environment.

For X-ray CT examples, the projection and back projection operator
subroutines may be called from MATLAB. To prepare MATLAB recognizable, go
to `imgRecSrc/prj` and execute `make mCPUPrj` for CPU implementation and
`make mGPUPrj` for GPU implementation of the X-ray CT projection operators.


### References
