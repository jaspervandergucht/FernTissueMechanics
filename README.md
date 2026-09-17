# FernTissueMechanics
Matlab code for Finite element analysis of fern tissue 

###################################################
Matlab code used for finite element calculations in:
S. Woudenberg, A.R.G. Plackett, Z. Hao, H. Suzuki, L. Alonso Baez, C. Borassi, T. Hamann, M. Ueda,  J.A. Langdale, J. Sprakel, J. van der Gucht, D. Weijers, "Transgenerational polarity axis inheritance during Ceratopteris embryogenesis", 2026 
by Jasper van der Gucht (Wageningen University)
###################################################

The code consists of two parts:
1. Analyze a bmp image of a plant tissue and extract the network of cell walls, and also create a triangular mesh in each cell [1].
2. Calculate stress distribution in the tissue for given stiffness gradient and turgor pressure using finite element methods

The matlab script main.m runs all stepos of the analysis for the provided .bmp image 

[1] Triangulation is done using Distmesh: P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB. SIAM Review, Volume 46 (2), pp. 329-345, 2004
