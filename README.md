# beta VDR
A matlab implementation of the methods presented in the paper *A Stable Finite Difference Method Based on Upward Continuation to Evaluate Vertical Derivatives of Potential Field Data*, published in Pure and Applied Geophysics

[https://doi.org/10.1007/s00024-022-03164-z]

**Authors:**  
* Saulo P. Oliveira [saulopo@ufpr.br]
* Luan Thanh Pham  [luanpt@hus.edu.vn]
         
## Instructions:

1. Make sure you have all the files provided in this package in the same folder:

* Example.m    : Main script  
* LoadXYZ.m    : Script to load the data
* taper2d.m    : Auxiliary function for 2D grid expansion   
* betaVDR.m    : Beta-VDR fourth-order finite-difference method 
* TN4th.m      : Fourth-order Tran-Nguyen method
* BDF.m        : First-order backward finite difference method
* dzFreq.m     : Vertical derivative in frequency domain
* input.dat    : Example of input data
* ColorMap.mat : colormap used in the manuscript
* license.txt  : License term
* readme.txt   : This file

2. If you prefer to use your own data, please save it as an ASCII file with 3 columns with the following data per column:

* Column 1: x coordinates 
* Column 2: y coordinates 
* Column 3: anomaly

If necessary, update the file name on Line 15 of file Example.m:

[x,y,Mo]=LoadXYZ('input.dat');

3. Run the program 'Example'

4. The program will generate, using the four methods above (betaVDR, TN4th, BDF, and dzFreq), the following three plots:

* Vertical derivative
* TDX (Cooper & Cowan, Comp & Geosc, 2006)
* ED (Pham, J Afr Eath Sci, 2021)
