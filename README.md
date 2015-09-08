## Image Reconstruction Source Code

### What does this package contain

#### 1. Accelerated Proximal-Gradient Algorithms (NPG and NPGs)

The examples about this methods are under `npgEx` folder with filenames
ended with `Ex`, where the code can reproduce the figures appear in our
paper.  The algorithm implementations are under folder `npg`, which use the
some utility functions under `utils` folder.

##### References

1. R. Gu and A. Dogandžić. (Feb. 2015). Reconstruction of nonnegative
   sparse signals using accelerated proximal-gradient algorithms. *arXiv*:
   [1502.02613](http://arxiv.org/abs/1502.02613).

<!---
R. Gu and A. Dogandžić, “Nesterov’s Proximal-Gradient Algorithms for Reconstructing Nonnegative Signals with Sparse Transform Coefficients,” 2014.
--->


####2. Beam Hardening Correction Algorithms

All examples and data for our blind beam hardening correction method are
under folder `bhcEx` with filenames ended with `Ex`.  The figures in our
paper can be reproduced by first run `*Ex`, e.g., `yangEx`, followed by
`*Ex('plot')`, e.g., `yangEx('plot')`.  Algorithm implementations are under
`bhc`.

##### References

1. R. Gu and A. Dogandžić, Beam hardening correction via mass attenuation
   discretization, in *Proc. IEEE Int. Conf. Acoust., Speech, Signal
   Process.*, Vancouver, Canada, May 2013, pp. 1085–1089.

1. R. Gu and A. Dogandžić, Polychromatic sparse image reconstruction and
   mass attenuation spectrum estimation via B-spline basis function
   expansion, in *Rev. Prog. Quant. Nondestr. Eval. Vol. 34 1650*, edited
   by D. E. Chimenti and L. J. Bond, AIP Conf. Proc. (2015), pp. 1707–1716.


### How to Install

To install this package, first download the repository by running

    git clone https://github.com/isucsp/imgRecSrc.git

after downloading, from MATLAB change your current folder to `imgRecSrc/`
and first execute `setupPath.m` to add necessary paths to the environment.

For X-ray CT examples, the projection and back projection operator
subroutines may be called from MATLAB.  Since they are written in `c`
language, to prepare MATLAB recognizable `MEX` files, go to `imgRecSrc/prj`
and compile the necessary files.  Instructions on compiling the code are
provided for both `UNIX` and `Windows`:

#### For `UNIX`

require: gcc, cuda toolkit (optional) and GPU (optional)

Execute `make cpu` to compile all cpu implementations.  If you have GPU
equipped, run `make gpu` to compile GPU implementation of the X-ray CT
projection operators.  The matlab code will automatically choose to run on
GPU if equipped.

If errors are reported while compiling the `*.c`/`*.cu` files under
`imgRecSrc/prj`, please edit the first few lines in
`imgRecSrc/prj/Makefile` to make sure the path for your `CUDA` installation
is correct.

#### For `Windows`

require: Visual Studio, cuda toolkit (optional) and GPU (optional)

Execute `imgRecSrc/setupPath.m` can automatically compile all needed files
and add paths.

If there is a GPU equipped in your PC, follow the following steps:

* Open the `VS Native Tools Command Prompt` via `Start -> Microsoft Visual
Studio -> Visual Studio Tools`;

* Use `cd` command to change directory to your `imgRecSrc/prj`;

* Run `nvcc -c gpuPrj.cu` to generate the `obj` file;

* Execute `imgRecSrc/setupPath.m`.

