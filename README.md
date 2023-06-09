[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8020304.svg)](https://doi.org/10.5281/zenodo.8020304)

# Reproduction Code

This repository contains the scripts to reproduce the numerical examples presented in the paper *Unfitted Trefftz discontinuous Galerkin methods for elliptic boundary value problems* by F. Heimann, C. Lehrenfeld, P. Stocker, and H. von Wahl.

### Installation
To run the python scripts locally, a compatible combination of `Netgen/NGSolve`, `ngsxfem` and `NGSTrefftz` are required. These can be installed by building from sources or the provided pip wheels. For detailed installation instructions, we refer to the installation guidelines of [NGSolve](https://docu.ngsolve.org/latest/install/install_sources.html), [ngsxfem](https://github.com/ngsxfem/ngsxfem/blob/release/INSTALLATION.md) and [NGSTrefftz](https://paulst.github.io/NGSTrefftz/readme.html#installing-the-package). Our numerical results are realised using the following versions:

| Package | git commit
|-|-|
| NGSolve | `42eddb00be6046078244b67a5a54c2bd3e598ad4`
| ngsxfem | `e3bcd370502ed24d8fc39a57cdbb1895276f6723`
| NGSTrefftz | `98289abd821ce1cd9a283bdb37cd0e3b652520d1`

To use the pip installer, please ensure that an appropriate pre-release versions of each package is installed.

### Content
| filename | description | 
|-|-|
| [`ficdom_dgT.py`](ficdom_dgT.py) | The python script to compute the unfitted Poisson problem using the unfitted DG/Trefftz methods, with or without geometry approximation and global or patch-wise ghost-penalties. |
| [`conv_study.py`](conv_study.py) | Python script to run a method for a given order over a series of meshes. |
| [`run_conv2d.sh`](run_conv2d.sh) | Shell script to run the two dimensional convergence studies. |
| [`run_conv3d.sh`](run_conv3d.sh) | Shell script to run the three dimensional convergence studies. |
| [`timings.py`](timings.py) | Python script to time the different methods. |
| [`run_timings.sh`](run_timings.sh) | Shell script to run the two and three dimensional timings. |
| `out/` | Directory containing the results realised with the given implementation. |
| [`ficdom_conv_diff_dgT.py`](ficdom_conv_diff_dgT.py) | Python script to compute the convection-diffusion problem with unfitted Trefftz/DG method. |

