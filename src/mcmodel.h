#ifndef __MCMODEL_H__
#define __MCMODEL_H__


/* Likelihood functions */
MrBFlt    LnProbAllocation (int *allocationVector, int numChars, MrBFlt alphaDir);
MrBFlt    LnProbEmission(int *latentPattern, int numCharsInCluster, int numMissing);
MrBFlt    LnProbLatentMatrix (int *allocationVector, int *latentMatrix, int numLatCols, int numChars, int compMatrixStart);

/* correlation model utility functions */
MrBLFlt   ExponentBySquaring(MrBLFlt base, MrBLFlt exp);
MrBLFlt   SmartExponentiation(MrBLFlt base, MrBLFlt exp);
int       GetNumLatentPatterns(int n, int d, int t);
int       *FlipCharacterPattern(int *pattern, int numCharsInCluster);
int       CheckCharacterPatternCompatibility(int *pattern, int *patternToCheck, int numCharsInCluster);
int       *GetClusterData(int *allocationVector, int cluster, int numCharsInCluster, int numChars, int compMatrixStart);
int       GetNumMissingChars(int numChars, int *pattern);
MrBFlt    *NormalizeEmissionProbabilities(MrBFlt *emissionProbabilities);
int       *ConvertLatentStates(int *dataSubset, int *origLatentPattern, int numCharsInCluster, int endStateIndex, MrBFlt rho, RandLong *seed, MrBFlt *moveProb, int forceOriginal);
int       *CountLatentResolutions(int *dataSubset, int *origLatentPattern, int numCharsInCluster, int endStateIndex, RandLong *seed);
int       *RescaleAllocationVector(int *allocationVector, int numChars, int newTable, int oldTable);
int       *DrawNewLatentMatrix(int *oldAllocationVector, int *newAllocationVector, int numChars, int compMatrixStart, int srcTable, int dstTable, int *oldLatentMatrix, MrBFlt rho, RandLong *seed, MrBFlt *moveProbSrc, MrBFlt *moveProbDst);
int       PrintLatentMatrixToScreen(int numChars, int *latentMatrix);
int       *DrawLatentPattern(int *data, int *origLatentStates, int numCharsInCluster, int numMissing, MrBFlt rho, RandLong *seed, MrBFlt *moveProb);
MrBFlt    ForceLatentPattern(int *data, int *origLatentStates, int numCharsInCluster, int numMissing, MrBFlt rho, RandLong *seed);
MrBFlt    ForceLatentPatternWithChangedData(int *oldAllocationVector, int *newAllocationVector, int *oldLatentMatrix, int *updatedLatentMatrix, int numChars, int compMatrixStart, int oldTable, int newTable, MrBFlt rho, RandLong *seed);


#endif  /* __MCMODEL_H__ */
