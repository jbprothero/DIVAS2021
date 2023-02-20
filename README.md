# DIVAS2021
MATLAB Code for Data Integration Via Analysis of Subspaces

Main function: DJIVEMainJP(), takes in the datablocks as a cell array, a classic Steve Marron paramstruct, and an optional third argument for "true signal" for diagnostic use.
Executes three subroutines:

  DJIVESignalExtractJP(), Finds signal ranks & perturbation angles for each data block. This function will attempt to initialize a parallel pool.
  
  DJIVEJointStrucEstimateJPLoadInfo(), Finds partially shared structure between data blocks. Runs an optimization problem implemented in CVX, which must be installed on your machine: http://cvxr.com/cvx/
  
  DJIVEReconstructMJ(), Solves for corresponding loadings structure to the scores structure found in the previous function. Packages relevant info into the returned "outstruct"
  
  outstruct components:
  
    jointBasisMap: Map data object with string integer keys for each block collection with found joint structure. Each container has an orthonormal basis matrix for the shared subspace in score space between the blocks for that container. Each array is n x K, where K is that collection's shared subspace rank.
    
    matLoadings: Cell array of map data objects containing the orthonormal basis matrices for "shared" subspaces in each block's loadings space. Each array is d_i x K, where d_i is the number of traits in block i and K is that collection's shared subspace rank.
    
    matBlocks: Cell array of map data objects containing signal modes of variation for each included block in each block collection. Each array is d_i x n, where d_i is the number of traits in block i.
    
    keyIdxMap: Map data object that converts between integer strings and block collection indices.
    rBars: Vector of the filtered rank of each data block
    
    phiBars: Vector of the scores-space perturbation angle of each data block
   
    psiBars: Vector of the loadings-space perturbation angle of each data block
    
    VBars: Cell array of orthonormal basis matrices for the estimated signal scores subspaces for each data block.
    
    UBars: Cell array of orthonormal basis matrices for the estimated signal loadings subspaces for each data block.
    
    VVHatCacheBars: Cell array of cell arrays of cached matrices for estimating diagnostic angle upper bounds in scores space.
    
    UUHatCacheBars: Cell array of cell arrays of cached matrices for estimating diagnostic angle upper bounds in loadings space.
  
Diagnostic function: DJIVEAngleDiagnosticJP(), takes in the datablocks and block names as cell arrays, the output from DJIVEMainJP(), an argument for random seed initialization, and a string for titling resulting figures. Creates three diagnostic images for angles between joint structure and estimated signal spaces. This function will attempt to initialize a parallel pool.

Please see the Google Drive folder linked below the input data for the demo script along with the output of DJIVEMainJP() on my machine for comparison: https://drive.google.com/drive/folders/1ToBIzVTfR2Hm53lCVneR0iGGTINR9r1w?usp=share_link
