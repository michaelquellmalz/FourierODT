FourierODT
=========================

Overview
--------
***FourierODT*** is a Matlab toolbox for optical diffraction tomography (ODT).
It implements reconstruction algorithms and phase retrieval methods 
for optical diffraction tomography.

__Reconstruction Methods__
* Backpropagation algorithm
* Inverse NDFT (nonuniform discrete Fourier transform) with the conjugate gradients method *(CGNE)*
* TV regularization with a forward-backward primal-dual algorithm *(FBPD)*

__Phase Retrieval Methods__
* Error reduction algorithm *(ER)*
* Hybrid input-output algorithm *(HIO)*
* Free space propagation approach by Maleki & Devaney *(MD)*

Reference
---------
When you are using this code, please cite the paper

Robert Beinert, Michael Quellmalz:
''Total variation-based phase retrieval for diffraction tomography''.
[ArXiv Preprint 2201.11579](https://arxiv.org/abs/2201.11579), 2022.

This paper also explains the algorithms in more detail.

Installation
------------
This software requires an installation of the NFFT3 Matlab interface 
in the Matlab path (by default add `../nfft`). 
This folder should contain nfft.m and nfftmex.mexw64.
The NFFT3 Matlab inteface is availabe at 
https://tu-chemnitz.de/~potts/nfft/download.php.
Then run `startup.m` in Matlab and you can execute the srcipts from `examples`.

Directory structure
-------------------

File/Folder        | Purpose
------------------:| ------------------------------------------------------
data_sets (dir)    | Contains data sets that matlab creates (initially empty)
examples (dir)     | Code for creating the figures in the paper
src (dir) 	       | Source code
COPYING            | License information
README.md          | This file
startup.m          | Startup Matlab script

Feedback
--------
Your comments are welcome! This is the first version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions in Guthub Issues.
Alternatively, you might contact
[Robert Beinert](mailto:beinert@math.tu-berlin.de)
or
[Michael Quellmalz](mailto:quellmalz@math.tu-berlin.de).

Legal Information & Credits
---------------------------
Copyright (c) 2022 Robert Beinert and Michael Quellmalz

This software was written by [Robert Beinert](https://www.math.tu-berlin.de/fachgebiete_ag_modnumdiff/angewandte_mathematik/v_menue/team/dr_robert_beinert/v_menue/dr_robert_beinert/) and [Michael Quellmalz](https://page.math.tu-berlin.de/~quellm/).
It was developed at the Institute of Mathematics, TU Berlin.
The second mentioned author acknowledges support by the German Research Foundation within the [SFB Tomography Across the Scales]( https://tomography.univie.ac.at/).

FourierODT is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. If not stated otherwise, this applies to all files contained in this
package and its sub-directories.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

