## Projected Nesterov's Proximal-Gradient Algorithm

The examples about this methods are under `npgEx` folder with filenames
ended with `Ex`, where the code can reproduce the figures appear in our
papers and reports.  The algorithm implementations are under folder `npg`,
which use the some utility functions under `utils` folder.

### How to Install

To install this package, first download the repository by running

    git clone https://github.com/isucsp/pnpg.git

after downloading, from MATLAB change your current folder to `pnpg/`
and execute `setupPath.m` to add necessary paths to the environment.

For `Windows`, you may need to have Visual Studio or other C/C++ compilers 
installed to compile some C code while calling `setupPath.m`.

For `UNIX`, you may need to have gcc installed. 

For the 3rd party softwares that are used, please refer to 
`getOthersCode.sh` in how to get them. 

### Note

The comments in some of `*.m` files may contain greek letters, which
are `UTF-8` encoded.  Please open in an appropriately configured text
editor.

### References

1. R. Gu and A. Dogandžić, (May. 2016). Projected Nesterov's
   Proximal-Gradient Algorithm for Sparse Signal Reconstruction with a
   Convex Constraint. *arXiv*: [1502.02613](http://arxiv.org/abs/1502.02613) \[stat.CO\].
   
1. R. Gu and A. Dogandžić, “Projected Nesterov’s proximal-gradient signal
   recovery from compressive Poisson measurements”, in *Proc. Asilomar Conf.
   Signals, Syst. Comput.*, Pacific Grove, CA, Nov. 2015, pp. 1490–1495.
   [\[DOI\]](http://dx.doi.org/10.1109/ACSSC.2015.7421393)
   [\[PDF\]](http://isucsp.github.io/pnpg/pdf/asilomar2015.pdf)

<!---
R. Gu and A. Dogandžić, “Nesterov’s Proximal-Gradient Algorithms for Reconstructing Nonnegative Signals with Sparse Transform Coefficients,” 2014.
--->


