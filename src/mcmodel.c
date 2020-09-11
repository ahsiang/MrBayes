/*
 *  MrBayes 3
 *
 *  (c) 2002-2013
 *
 *  John P. Huelsenbeck
 *  Dept. Integrative Biology
 *  University of California, Berkeley
 *  Berkeley, CA 94720-3140
 *  johnh@berkeley.edu
 *
 *  Fredrik Ronquist
 *  Swedish Museum of Natural History
 *  Box 50007
 *  SE-10405 Stockholm, SWEDEN
 *  fredrik.ronquist@nrm.se
 *
 *  With important contributions by
 *
 *  Paul van der Mark (paulvdm@sc.fsu.edu)
 *  Maxim Teslenko (maxkth@gmail.com)
 *  Chi Zhang (zhangchicool@gmail.com)
 *
 *  and by many users (run 'acknowledgments' to see more info)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 */

 #include "bayes.h"
 #include "best.h"
 #include "command.h"
 #include "mcmc.h"
 #include "model.h"
 #include "utils.h"
 #include "likelihood.h"
 #include "mcmodel.h"
 #include "proposal.h"

 #if !defined(MAX)
 #define MAX(a,b)                            (((a) > (b)) ? (a) : (b))
 #endif

 /* local prototypes */
 MrBLFlt Combination (int n, int k);


 /*---------------------------------------------------------------------------------
 |
 |   LnProbEmission
 |
 |   Calculates log emission probability for a latent pattern.
 |   Missing character compliant.
 |
 ---------------------------------------------------------------------------------*/
 MrBFlt LnProbEmission(int *latentPattern, int numCharsInCluster, int numMissing)
 {
     int            i, n, m=0, q, p;
     MrBFlt         lnProbability;

     /* Calculate total emission probability */
     n = numCharsInCluster;
     q = numMissing;

     /* Count number of i-states in latent pattern */
     for (i=0; i<numLocalTaxa; i++)
         if (latentPattern[i] == INTSTATE)
             m++;

     /* Get mn-q power term */
     p = m * n - q;

     /* Get emission probability */
     lnProbability = p * (log(1) - log(2)); // i.e., (1/2)^(mn-q)

     return ( lnProbability );
 }


 /*-----------------------------------------------------------------
 |
 |   LnProbLatentMatrix: Calculate likelihood of latent matrix
 |
 -----------------------------------------------------------------*/
 MrBFlt LnProbLatentMatrix (int *allocationVector, int *latentMatrix, int numClusters, int numChars, int compMatrixStart)
 {
     int         i, j, idx=0, highest=0, clusterCols[numClusters], currColumn[numLocalTaxa],
                 numCharsPerCluster[numClusters], numMissingPerCluster[numClusters], numDataValues,
                 *data;
     MrBFlt      lnTotalProb=0.0;

     /* Determine which columns to calculate likelihood for */
     for (i=0; i<numChars; i++)
         if (allocationVector[i] == highest)
             {
             clusterCols[idx++] = i;
             highest++;
             }

     /* Get number of characters per cluster */
     for (i=0; i<numClusters; i++)
         numCharsPerCluster[i] = 0; // Initialize all entries to 0
     for (i=0; i<numClusters; i++)
         {
         for (j=0; j<numChars; j++)
             if (allocationVector[j] == i)
                 numCharsPerCluster[i]++;
         }

     /* Get number of missing characters per cluster */
     for (i=0; i<numClusters; i++)
         numMissingPerCluster[i] = 0;
     for (i=0; i<numClusters; i++)
         {
         // Grab relevant data
         data = GetClusterData(allocationVector, i, numCharsPerCluster[i], numChars, compMatrixStart);
         // Count number of missing characters
         numDataValues = numLocalTaxa * numCharsPerCluster[i];
         for (j=0; j<numDataValues; j++)
             if (data[j] == MISSING || data[j] == GAP)
                 numMissingPerCluster[i]++;
         // Free allocation
         free (data);
         data = NULL;
         }

     /* Loop through processes and calculate probability */
     for (i=0; i<numClusters; i++)
         {
         // Grab current column/cluster/process
         for (j=0; j<numLocalTaxa; j++)
             currColumn[j] = latentMatrix[pos(j, clusterCols[i], numChars)];

         // Calculate log likelihood of column and add to total
         lnTotalProb += LnProbEmission(currColumn, numCharsPerCluster[i], numMissingPerCluster[i]);
         }

     return lnTotalProb;
 }



/*---------------------------------------------------------------------------------
|
|   ExponentBySquaring
|
|   Carries out exponentiation using the exponentiation by squaring algorithm.
|
---------------------------------------------------------------------------------*/
MrBLFlt ExponentBySquaring(MrBLFlt base, MrBLFlt exp)
{
    if (exp < 0)
        return ( ExponentBySquaring(1.0/base, -exp) );
    else if ((int) exp == 0)
        return ( 1 );
    else if ((int) exp == 1)
        return ( base );
    else if ((int) exp % 2 != 0) // Odd exponent case
        return ( base * ExponentBySquaring(base*base, (exp-1.0)/2.0) );
    else // Even exponent case
        return ( ExponentBySquaring(base*base, exp/2.0) );
}


/*---------------------------------------------------------------------------------
|
|   SmartExponentiation
|
|   Wrapper function for exponentation which picks the best method depending on the
|   magnitude of the exponentiation required.
|
---------------------------------------------------------------------------------*/
MrBLFlt SmartExponentiation(MrBLFlt base, MrBLFlt exp)
{
    int       i, b, e;
    MrBLFlt   result=1.0;

    b = (int) base;
    e = (int) exp;

    if ((b == 2) && (e < 64)) // Use bitshifting for smaller powers of 2
        result = (MrBLFlt) (1ULL << (BitsLong) exp);
    else
        {
        if (e < 64) // Use naive looping for small values of exp
            {
            for (i=0; i<e; i++)
                result *= base;
            }
        else // Use exponentiation by squaring for large values of exp
            result = ExponentBySquaring(base, exp);
        }

    return (result);
}


/*---------------------------------------------------------------------------------
|
|   GetClusterData
|
|   Grabs the character data in a particular cluster.
|
---------------------------------------------------------------------------------*/
int *GetClusterData(int *allocationVector, int cluster, int numCharsInCluster, int numChars, int compMatrixStart)
{

    int         i, j, numDataValues, *data, idx=0;

    /* Allocate space for data */
    numDataValues = numLocalTaxa * numCharsInCluster;
    data = (int *) SafeMalloc ((size_t)numDataValues * sizeof(int));
    if (!data)
        printf("ERROR: Problem with allocation in GetClusterData\n");

    for (i=0; i<numChars; i++)
        {
        if (allocationVector[i] == cluster)
            {
            for (j=0; j<numLocalTaxa; j++)
                data[pos(j, idx, numCharsInCluster)] = compMatrix[pos(j, compMatrixStart+i, numCompressedChars)];
            idx++;
            }
        }

    return (data);
}


/*---------------------------------------------------------------------------------
|
|   FlipCharacterPattern
|
|   Takes a character pattern and finds its opposite, e.g., 000 -> 111.
|
---------------------------------------------------------------------------------*/
int *FlipCharacterPattern(int *pattern, int numCharsInCluster)
{
    int         i, *flippedPattern;

    /* Allocate space for flipped pattern */
    flippedPattern = (int *) SafeMalloc ((size_t)(numCharsInCluster) * sizeof(int));
    if (!flippedPattern)
        printf("ERROR: Problem with allocation in FlipCharacterPattern\n");

    for (i=0; i<numCharsInCluster; i++)
        {
        if ((pattern[i] == MISSING) || (pattern[i] == GAP))
            flippedPattern[i] = pattern[i];
        else
            {
            if (pattern[i] == 1)
                flippedPattern[i] = 2;
            else
                flippedPattern[i] = 1;
            }
        }

    return (flippedPattern);
}


/*---------------------------------------------------------------------------------
|
|   GetNumLatentPatterns
|
|   Calculates the number of latent patterns possible for a given number of
|   i-states (n) given the number of dimorphic (d; i.e., 0/i or i/1) and trimorphic
|   (t; i.e., 0/i/1) states.
|
---------------------------------------------------------------------------------*/
int GetNumLatentPatterns(int n, int d, int t)
{
    int         sum=0, term, i;


    for (i=0; i<n+1; i++)
        {
        term = MAX(1, 2*(t-(n-i)));
        sum += Combination(d,i) * Combination(t,n-i) * term;
        }

    return (sum);

}


/*---------------------------------------------------------------------------------
|
|   CheckCharacterPatternCompatibility
|
|   Determines if a character pattern is 0- or 1-compatible with another specified
|   character pattern.
|
|   Returns ENDSTATE if 0-compatible, OPPENDSTATE if 1-compatible, INTSTATE
|   if neither, and TRIPOLY if completely ambiguous.
|
---------------------------------------------------------------------------------*/
int CheckCharacterPatternCompatibility(int *pattern, int *patternToCheck, int numCharsInCluster)
{
    int         i, *oppositePattern, numAmbiguous=0, equal=TRUE, opposite=TRUE, result;


    /* Get opposite pattern of specified character pattern */
    oppositePattern = FlipCharacterPattern(pattern, numCharsInCluster);

    /* Loop through characters to check compatibility */
    for (i=0; i<numCharsInCluster; i++)
        {
        if ((pattern[i] == MISSING) || (pattern[i] == GAP) || (patternToCheck[i] == MISSING) || (patternToCheck[i] == GAP))
            {
            numAmbiguous++;
            continue;
            }
        else
            {
            if (patternToCheck[i] != pattern[i])
                equal = FALSE;
            if (patternToCheck[i] != oppositePattern[i])
                opposite = FALSE;
            }
        if ((equal == FALSE) && (opposite == FALSE))
            break;
        }

    /* Determine compatibility */
    if (numAmbiguous == numCharsInCluster)
        result = TRIPOLY;
    else
        {
        if (equal == TRUE)
            result = ENDSTATE;
        else
            {
            if (opposite == TRUE)
                result = OPPENDSTATE;
            else
                result = INTSTATE;
            }
        }

    // Free allocations
    free (oppositePattern);
    oppositePattern = NULL;

    return (result);

}


/*---------------------------------------------------------------------------------
|
|   CountLatentResolutions
|
|   Takes a subset of data that comprises a single cluster and determines
|   the number of latent resolutions for each possible number of i-states,
|   given an end state index. Returns a vector [r_0,r_1,...,r_k], where k
|   is the maximum number of i-states possible (= number of taxa), and r_k
|   is the number of possible latent resolutions that contain k i-states.
|
|   Uses the following codes for the possible latent states:
|       0       ENDSTATE
|       1       OPPENDSTATE
|       i       INTSTATE
|       0/i/1   TRIPOLY
|       0/i     ZEROPOLY
|       i/1     ONEPOLY
|
---------------------------------------------------------------------------------*/
int *CountLatentResolutions(int *dataSubset, int *origLatentPattern, int numCharsInCluster, int endStateIndex, RandLong *seed)
{
    int         i, j, *latentResolutionCounts, endStatePattern[numCharsInCluster],
                possibleLatentStates[numLocalTaxa], currentCharacterPattern[numLocalTaxa],
                charCompat, numStaticIntStates=0, numDiPolys=0, numTriPolys=0;


    /* Allocate space for vector of latent resolution counts */
    latentResolutionCounts = (int *) SafeMalloc ((size_t)(numLocalTaxa+1) * sizeof(int));
    if (!latentResolutionCounts)
        printf("ERROR: Problem with allocation in CountLatentResolutions\n");

    /* Get end state pattern */
    for (i=0; i<numCharsInCluster; i++)
        endStatePattern[i] = dataSubset[pos(endStateIndex, i, numCharsInCluster)];

    /* Loop through character patterns to determine possible latent states */
    for (i=0; i<numLocalTaxa; i++)
        {
        // Specified end state is trivial
        if (i == endStateIndex)
            possibleLatentStates[i] = ENDSTATE;
        else
            {
            // First get the current character pattern
            for (j=0; j<numCharsInCluster; j++)
                currentCharacterPattern[i] = dataSubset[pos(i, j, numCharsInCluster)];

            // Determine if this character pattern is equal to the end state or the opposite end state
            charCompat = CheckCharacterPatternCompatibility(endStatePattern, currentCharacterPattern, numCharsInCluster);

            // Set code for possible latent states based on pattern compatibility
            // Also count how many of each type of state there is
            switch(charCompat)
                {
                case ENDSTATE:
                    possibleLatentStates[i] = ZEROPOLY;
                    numDiPolys++;
                case OPPENDSTATE:
                    possibleLatentStates[i] = ONEPOLY;
                    numDiPolys++;
                case INTSTATE:
                    possibleLatentStates[i] = INTSTATE;
                    numStaticIntStates++;
                case TRIPOLY:
                    possibleLatentStates[i] = TRIPOLY;
                    numTriPolys++;
                }
            }
        }

    /* Loop through possible latent states to determine number of latent patterns per i-state */
    for (i=0; i<numLocalTaxa+1; i++)
        {
        // The last entry is for the single all-i-state
        if (i == numLocalTaxa)
            latentResolutionCounts[i] = 1;
        else
            {
            // The number of resolutions with n i-states is 0 for all n < numStaticIntStates
            if (i < numStaticIntStates)
                latentResolutionCounts[i] = 0;
            // In all other cases, we need to count how many resolutions there are
            else
                {
                if (i == numStaticIntStates)
                    latentResolutionCounts[i] = 1;
                // Remaining cases are when num i-states > numStaticIntStates
                else
                    latentResolutionCounts[i] = GetNumLatentPatterns(i-numStaticIntStates, numDiPolys, numTriPolys);
                }
            }
        }

    return (latentResolutionCounts);

}


/*---------------------------------------------------------------------------------
|
|   NormalizeEmissionProbabilities
|
|   Takes a set of log emission probabilities and returns normalized probabilities.
|
---------------------------------------------------------------------------------*/
MrBFlt *NormalizeEmissionProbabilities(MrBFlt *emissionProbabilities)
{

    int         i;
    MrBFlt     *normProbs, maxLogProb, differences[numLocalTaxa+1], threshProbs[numLocalTaxa+1],
                threshold, epsilon=1.0e-16, probSum=0.0;

    /* Allocate space for normalized probabilities */
    normProbs = (MrBFlt *) SafeMalloc ((size_t)(numLocalTaxa+1) * sizeof(MrBFlt));
    if (!normProbs)
        printf("ERROR: Problem with allocation in NormalizeEmissionProbabilities\n");

    // Get maximum log probability
    maxLogProb = emissionProbabilities[0];
    for (i=1; i<numLocalTaxa+1; i++)
        if (emissionProbabilities[i] > maxLogProb)
            maxLogProb = emissionProbabilities[i];

    // Subtract maximum from all log probs
    for (i=0; i<numLocalTaxa+1; i++)
        differences[i] = emissionProbabilities[i] - maxLogProb;

    // Compare differences to a threshold and only keep those that exceed it
    threshold = log(epsilon) - log(numLocalTaxa+1);
    for (i=0; i<numLocalTaxa+1; i++)
        {
        if (differences[i] < threshold)
            threshProbs[i] = 0;
        else
            threshProbs[i] = exp(differences[i]);
        }

    // Normalize remaining probabilities
    for (i=0; i<numLocalTaxa+1; i++)
        probSum += threshProbs[i];
    for (i=0; i<numLocalTaxa+1; i++)
        normProbs[i] = threshProbs[i] / probSum;

    return (normProbs);
}


/*---------------------------------------------------------------------------------
|
|   ConvertLatentStates
|
|   Takes a subset of data that comprises a single cluster and determines
|   the possible latent pattern that would arise, given an end state index.
|   If end state index is negative, then the all-intermediate-state case is
|   returned. Missing character compliant. Used for latent move.
|
|   If forceOriginal == YES, then this function will return the original latent
|   pattern and the corresponding move probability for that pattern.
|
---------------------------------------------------------------------------------*/
int *ConvertLatentStates(int *dataSubset, int *origLatentPattern, int numCharsInCluster, int endStateIndex, MrBFlt rho, RandLong *seed, MrBFlt *moveProb, int forceOriginal)
{
    int         i, j, *latentResolution, endStatePattern[numCharsInCluster],
                numMissingEndState, currentCharacterPattern[numCharsInCluster],
                compatibleLatentState, selectedLatentState;
    MrBFlt      emProb, iProb, jProb, pis[3], rand;


    /* Allocate space for final latent resolution */
    latentResolution = (int *) SafeMalloc ((size_t)numLocalTaxa * sizeof(int));
    if (!latentResolution)
        printf("ERROR: Problem with allocation in ConvertLatentStates\n");

    /* Set end state (or all i-state pattern) */
    if (endStateIndex == -1)
        {
        // End state pattern is unknown
        for (i=0; i<numCharsInCluster; i++)
            endStatePattern[i] = MISSING;
        }
    else
        {
        latentResolution[endStateIndex] = ENDSTATE;
        /* Get new end state pattern */
        for (i=0; i<numCharsInCluster; i++)
            endStatePattern[i] = dataSubset[pos(endStateIndex, i, numCharsInCluster)];
        }

    // if (numCharsInCluster > 1)
    // {
    //     printf("end state index: %d\n", endStateIndex);
    //     printf("end state pattern: ");
    //     for (i=0; i<numCharsInCluster; i++)
    //         printf("%c ",WhichStand(endStatePattern[i]));
    //     printf("\n");
    // }

    /* We need to determine the probability of end state/opposite end state patterns being emitted by an i-state. */
    // First get the number of missing states in the end state pattern.
    numMissingEndState = GetNumMissingChars(numCharsInCluster, endStatePattern);

    // Now get the emission probability of an i-state emitting this new end state pattern
    MrBFlt power = numCharsInCluster - numMissingEndState;
    emProb = SmartExponentiation(0.5, power);

    // if (numCharsInCluster > 1)
    // {
    //     printf("num missing end state: %d\n",numMissingEndState);
    //     printf("em prob: %f\n",emProb);
    // }

    // Now get the overall probability of an i-state emitting the end or opposite end state
    // (i.e., emProb * stationary state frequency of the i-state). First compute pis:
    pis[0] = pis[2] = 1.0 / (2.0 + rho);   // alpha
    pis[1] = 1.0 - (2.0 * pis[0]);         // beta
    // Finally get the overall probability
    iProb = log(emProb) + log(pis[1]);

    // Also get the probability of the end or opposite end state emitting the pattern
    jProb = log(1.0 - (emProb * pis[1]));

    // if (numCharsInCluster > 1)
    // {
    //     printf("iProb: %f\n",iProb);
    //     printf("jProb: %f\n",jProb);
    // }

    /* All intermediate case */
    if (endStateIndex < 0)
        {
        for (i=0; i<numLocalTaxa; i++)
            {
            // Set latent resolution to i-state
            latentResolution[i] = INTSTATE;

            // In the all i-state case, we don't know what the end state pattern is
            // So the probability of the underlying pattern matching any of the states
            // is equal to 1/3
            *moveProb += log(1.0/3.0);
            }
        }
    /* All other cases */
    else
        {
        /* Loop through original latent states to determine new latent states */
        for (i=0; i<numLocalTaxa; i++)
            {
            if (i == endStateIndex)
                continue;
            else
                {
                // If we don't force the original latent states, then we draw new ones probabilistically

                // Get current character pattern
                for (j=0; j<numCharsInCluster; j++)
                    currentCharacterPattern[j] = dataSubset[pos(i, j, numCharsInCluster)];

                // if (numCharsInCluster > 2)
                // {
                //     printf("current char pattern: ");
                //     for (j=0; j<numCharsInCluster; j++)
                //         printf("%c ",WhichStand(currentCharacterPattern[j]));
                //     printf("\n");
                // }

                // Get possible latent state(s)
                compatibleLatentState = CheckCharacterPatternCompatibility(endStatePattern, currentCharacterPattern, numCharsInCluster);
                //
                // if (numCharsInCluster > 2)
                // {
                //     printf("compat: %d\n",compatibleLatentState);
                // }

                if (forceOriginal == YES)
                    selectedLatentState = origLatentPattern[i];

                // Select latent states as appropriate
                if (compatibleLatentState == ENDSTATE)
                    {
                    // The underlying latent state could be 0 or i
                    // Determine if we pick the i-state
                    if (forceOriginal == NO)
                        {
                        rand = log(RandomNumber(seed));
                        if (rand <= iProb)
                            selectedLatentState = INTSTATE;
                        else
                            selectedLatentState = ENDSTATE;
                        }

                    // Now deal with the moveprob
                    if (selectedLatentState == INTSTATE)
                        *moveProb += iProb;
                    else
                        *moveProb += jProb;
                    }
                else if (compatibleLatentState == OPPENDSTATE)
                    {
                    // The underlying latent state could be 1 or i
                    // Determine if we pick the i-state
                    if (forceOriginal == NO)
                        {
                        rand = log(RandomNumber(seed));
                        if (rand <= iProb)
                            selectedLatentState = INTSTATE;
                        else
                            selectedLatentState = OPPENDSTATE;
                        }

                    // Now deal with the moveprob
                    if (selectedLatentState == INTSTATE)
                        *moveProb += iProb;
                    else
                        *moveProb += jProb;
                    }
                else if (compatibleLatentState == INTSTATE)
                    {
                    // The underlying latent state can only be i
                    selectedLatentState = INTSTATE;
                    }
                else
                    {
                    // The underlying state can be 0, i, or 1
                    // Should this be equal probabilities or equal to stationary state frequencies?
                    if (forceOriginal == NO)
                        {
                        rand = RandomNumber(seed);
                        if (rand <= 1.0/3.0)
                            selectedLatentState = ENDSTATE;
                        else if ((1.0/3.0 < rand) && (rand <= 2.0/3.0))
                            selectedLatentState = INTSTATE;
                        else
                            selectedLatentState = OPPENDSTATE;
                        }

                    *moveProb += log(1.0/3.0);
                    }

                latentResolution[i] = selectedLatentState;
                }
            }
        }

    return (latentResolution);
}


/*---------------------------------------------------------------------------------
|
|   RescaleAllocationVector
|
|   Rescales allocation vector to follow growth function.
|
---------------------------------------------------------------------------------*/
int *RescaleAllocationVector(int *allocationVector, int numChar, int newTable, int oldTable)
{
    int         i, j, barrier=-1, toSwitch;

    int *rescaled = SafeMalloc(numChar * sizeof(int*));
    if (!rescaled)
        printf("ERROR: Problem with allocation in RescaleAllocationVector\n");

    for (i=0; i<numChar; i++)
        rescaled[i] = allocationVector[i];

    if (newTable != oldTable)
        {
        for (i=0; i<numChar; i++)
            {
            if (rescaled[i] == barrier + 1)
                barrier++;
            else
                {
                if (rescaled[i] > barrier + 1)
                    {
                    toSwitch = rescaled[i];
                    for (j=i; j<numChar; j++)
                        {
                        if (rescaled[j] == toSwitch)
                            rescaled[j] = barrier + 1;
                        else if (rescaled[j] >= barrier + 1)
                            rescaled[j]++;
                        }
                    barrier++;
                    }
                }
            }
        }
    return (rescaled);
}


/*---------------------------------------------------------------------------------
|
|   DrawNewLatentMatrix
|
|   Draws new latent matrix that match data pattern for new clusters.
|
---------------------------------------------------------------------------------*/
int *DrawNewLatentMatrix(int *oldAllocationVector, int *newAllocationVector, int numChars, int compMatrixStart, int srcTable, int dstTable, int *oldLatentMatrix, MrBFlt rho, RandLong *seed, MrBFlt *moveProbSrc, MrBFlt *moveProbDst)
{
    int         i, j, numCharsInClusterDst=0, *dataDst,
                numMissingDst, *newLatentStatesDst,
                *finalLatentMatrix, numValues,
                numCharsInClusterSrc=0, drawSrc=FALSE, *dataSrc, numMissingSrc,
                *newLatentStatesSrc;

    /* Source cluster */
    /* Get number of characters in source cluster */
    for (i=0; i<numChars; i++)
        if (newAllocationVector[i] == srcTable)
            numCharsInClusterSrc++;

    if (numCharsInClusterSrc > 0) // There is a chance that the source cluster has disappeared
    {
        drawSrc = TRUE;

        /* Get data from characters that belong to the source cluster */
        dataSrc = GetClusterData(newAllocationVector, srcTable, numCharsInClusterSrc, numChars, compMatrixStart);

        /* Get number of missing characters */
        numMissingSrc = GetNumMissingChars(numCharsInClusterSrc * numLocalTaxa, dataSrc);

        /* Set new latent states */
        newLatentStatesSrc = DrawLatentPattern(dataSrc, NULL, numCharsInClusterSrc, numMissingSrc, rho, seed, moveProbSrc);

        // printf("src data:\n");
        // for (i=0; i<numLocalTaxa; i++)
        // {
        //     for (int j=0; j<numCharsInClusterSrc; j++)
        //         printf("%c ", WhichStand(dataSrc[pos(i,j,numCharsInClusterSrc)]));
        //     printf("\n");
        // }
        // printf("\n");
        //
        // printf("new latent states src: ");
        // for (i=0; i<numLocalTaxa; i++)
        //     printf("%c ",WhichLatent(newLatentStatesSrc[i]));
        // printf("\n");
        //
        // printf("move prob src: %f\n",*moveProbSrc);
        //
        // free (dataSrc);
        // dataSrc = NULL;
    }
    else
    {
        *moveProbSrc = 1.0;
    }

    /* Destination cluster */
    /* Get number of characters in destination cluster*/
    for (i=0; i<numChars; i++)
        if (newAllocationVector[i] == dstTable)
            numCharsInClusterDst++;

    /* Get data from characters that belong to the newTable cluster */
    dataDst = GetClusterData(newAllocationVector, dstTable, numCharsInClusterDst, numChars, compMatrixStart);

    // printf("data dst: \n");
    // for (i=0; i<numLocalTaxa; i++)
    // {
    //     for (j=0; j<numCharsInClusterDst; j++)
    //         printf("%c ",WhichStand(dataDst[pos(i,j,numCharsInClusterDst)]));
    //     printf("\n");
    // }
    // printf("\n");


    /* Get number of missing characters */
    numMissingDst = GetNumMissingChars(numCharsInClusterDst*numLocalTaxa, dataDst);


    /* Set new latent states */
    newLatentStatesDst = DrawLatentPattern(dataDst, NULL, numCharsInClusterDst, numMissingDst, rho, seed, moveProbDst);

    // printf("new latent states dst: ");
    // for (i=0; i<numLocalTaxa; i++)
    //     printf("%c ",WhichLatent(newLatentStatesDst[i]));
    // printf("\n");
    //
    // printf("move prob dst: %f\n",*moveProbDst);


    /* Initialize data structure to hold new latent matrix */
    numValues = numChars * numLocalTaxa;
    finalLatentMatrix = (int *) SafeMalloc ((size_t)numValues * sizeof(MrBFlt));
    if (!finalLatentMatrix)
        printf("ERROR: Problem with allocation in UpdateLatentPatterns\n");

    /* Replace appropriate columns with new latent state resolution */
    for (i=0; i<numLocalTaxa; i++)
        {
        for (j=0; j<numChars; j++)
            {
            if ((newAllocationVector[j] == srcTable) && (drawSrc == TRUE))
                finalLatentMatrix[pos(i, j, numChars)] = newLatentStatesSrc[i];
            if (newAllocationVector[j] == dstTable)
                finalLatentMatrix[pos(i, j, numChars)] = newLatentStatesDst[i];
            else
                finalLatentMatrix[pos(i, j, numChars)] = oldLatentMatrix[pos(i, j, numChars)];
            }
        }

    /* Free allocations */
    if (drawSrc == TRUE)
    {
        free (newLatentStatesSrc);
        newLatentStatesSrc = NULL;
    }

    free (dataDst);
    dataDst = NULL;

    free (newLatentStatesDst);
    newLatentStatesDst = NULL;


    return (finalLatentMatrix);

}


/*---------------------------------------------------------------------------------
|
|  DrawLatentPattern
|
|  Takes a set of data associated with a given cluster and probabilistically
|  draws a latent resolution
|
---------------------------------------------------------------------------------*/
int *DrawLatentPattern(int *data, int *origLatentStates, int numCharsInCluster, int numMissing, MrBFlt rho, RandLong *seed, MrBFlt *moveProb)
{
    int         i, j, idx, *latentResolution, *newLatentPattern,
                newLatentVectorIdx=-1;
    MrBFlt      latentResolutions[(numLocalTaxa+1) * numLocalTaxa],
                *normalizedProbs, cumulativeProbs[numLocalTaxa+1], rand, currentMoveProb,
                allMoveProbs[numLocalTaxa+1], emProb;

    // Allocate space for final latent resolution
    newLatentPattern = (int *) SafeMalloc(numLocalTaxa * sizeof(int));
    if (!newLatentPattern)
        printf("ERROR: Problem with allocation in DrawLatentPattern\n");

    // Loop through all possible end states
    for (i=0; i<numLocalTaxa+1; i++)
        {
        // Reset currentMoveProb
        currentMoveProb = 0.0;

        // Get current index
        if (i == numLocalTaxa) // Last entry is the all-i-state case
            idx = -1;
        else // Current data pattern is the end state
            idx = i;

        // Get latent resolution
        latentResolution = ConvertLatentStates(data, origLatentStates, numCharsInCluster, idx, rho, seed, &currentMoveProb, NO);

        // Copy latent resolution to data structure
        for (j=0; j<numLocalTaxa; j++)
            latentResolutions[pos(j, i, numLocalTaxa+1)] = latentResolution[j];

        // Get emission probability of latent resolution and then save to data structure
        emProb = LnProbEmission(latentResolution, numCharsInCluster, numMissing);

        // Calculate total move prob for current resolution
        allMoveProbs[i] = emProb + currentMoveProb;

        // Free allocation
        free (latentResolution);
        latentResolution = NULL;
        }

    /* Normalize emission probabilities */
    normalizedProbs = NormalizeEmissionProbabilities(allMoveProbs);

    // printf("latent resolutions: \n");
    // for (i=0; i<numLocalTaxa; i++)
    // {
    //     for (j=0; j<numLocalTaxa+1; j++)
    //         printf("%c ",WhichLatent(latentResolutions[pos(i,j,numLocalTaxa+1)]));
    //     printf("\n");
    // }
    // printf("\n");

    /* Convert to cumulative probabilities */
    for (i=0; i<numLocalTaxa+1; i++)
        {
        if (i == 0)
            cumulativeProbs[i] = normalizedProbs[i];
        else
            cumulativeProbs[i] = normalizedProbs[i] + cumulativeProbs[i-1];
        }

    // printf("cumulative probs: ");
    // for (i=0; i<numLocalTaxa+1; i++)
    //     printf("%f ",cumulativeProbs[i]);
    // printf("\n");

    /* Randomly select new latent resolution */
    rand = RandomNumber(seed);
    for (i=0; i<numLocalTaxa+1; i++)
        {
        if (rand <= cumulativeProbs[i])
            {
            newLatentVectorIdx = i;
            break;
            }
        }

    // printf("new latent vector idx: %d\n",newLatentVectorIdx);

    /* Get new latent resolution pattern */
    for (i=0; i<numLocalTaxa; i++)
        newLatentPattern[i] = latentResolutions[pos(i, newLatentVectorIdx, numLocalTaxa+1)];

    // printf("new latent pattern: ");
    // for (i=0; i<numLocalTaxa; i++)
    //     printf("%c ",WhichLatent(newLatentPattern[i]));
    // printf("\n");

    /* Set move prob for the selected latent pattern */
    /* Note that this is NOT a log prob! */
    *moveProb = normalizedProbs[newLatentVectorIdx];
    //
    // printf("move prob (correct): %f\n",normalizedProbs[newLatentVectorIdx]);
    // printf("move prob: %f\n",*moveProb);
    // //
    // getchar();

    /* Free allocations */

    free (normalizedProbs);
    normalizedProbs = NULL;

    return newLatentPattern;

}



/*---------------------------------------------------------------------------------
|
|  ForceLatentPattern
|
|  Takes a set of data associated with a given cluster gets move probability
|  associated with the provided latent pattern. Used to get move prob for backwards
|  moves involving the latent matrix.
|
---------------------------------------------------------------------------------*/
MrBFlt ForceLatentPattern(int *data, int *origLatentStates, int numCharsInCluster, int numMissing, MrBFlt rho, RandLong *seed)
{
    int         i, idx, origEndStateIdx=-1, *latentResolution;
    MrBFlt      backwardsMoveProbs[numLocalTaxa+1], backwardsSecProb, backwardsEmProb,
                *backwardsNormProbs, forcedProb;

    // printf("----------------------\norig latent resolution: ");
    // for (int j=0; j<numLocalTaxa; j++)
    //     printf("%c ", WhichLatent(origLatentStates[j]));
    // printf("\n\n");
    //
    // printf("data:\n");
    // for (int j=0; j<numLocalTaxa; j++)
    // {
    //     for (int k=0; k<numCharsInCluster; k++)
    //         printf("%c ",WhichStand(data[pos(j,k,numCharsInCluster)]));
    //     printf("\n");
    // }
    //
    // printf("\n----------------------\n");

    // Get end state index in original latent pattern
    for (i=0; i<numLocalTaxa; i++)
        {
        if (origLatentStates[i] == ENDSTATE)
            {
            origEndStateIdx = i;
            break;
            }
        }
    if (origEndStateIdx == -1)
        origEndStateIdx = numLocalTaxa;

    // Get emission probabilities for every end state possibility
    for (i=0; i<numLocalTaxa+1; i++)
        {
        if (i == numLocalTaxa) // Last entry is the all-i-state case
            idx = -1;
        else // Current data pattern is the end state
            idx = i;

        // Reset move prob for current latent resolution
        backwardsSecProb = 0.0;

        // Get latent resolution
        if (i == origEndStateIdx)
            latentResolution = ConvertLatentStates(data, origLatentStates, numCharsInCluster, idx, rho, seed, &backwardsSecProb, YES);
        else
            latentResolution = ConvertLatentStates(data, origLatentStates, numCharsInCluster, idx, rho, seed, &backwardsSecProb, NO);

        backwardsEmProb = LnProbEmission(latentResolution, numCharsInCluster, numMissing);

        // printf("i: %d\n",i);
        // printf("latent resolution: ");
        // for (int j=0; j<numLocalTaxa; j++)
        //     printf("%c ", WhichLatent(latentResolution[j]));
        // printf("\n");
        // printf("backwards em prob: %f\n",backwardsEmProb);
        // printf("backwards sec prob: %f\n",backwardsSecProb);

        backwardsMoveProbs[i] = backwardsSecProb + backwardsEmProb;

        // Free allocations
        free (latentResolution);
        latentResolution = NULL;
        }


    /* Normalize emission probabilities */
    backwardsNormProbs = NormalizeEmissionProbabilities(backwardsMoveProbs);

    /* Get move probability for forced latent resolution */
    forcedProb = backwardsNormProbs[origEndStateIdx];

    // printf("forced prob: %f\n\n",forcedProb);
    // getchar();


    /* Free allocations */
    free (backwardsNormProbs);
    backwardsNormProbs = NULL;

    return forcedProb;
}



/*---------------------------------------------------------------------------------
|
|  ForceLatentPatternWithChangedData
|
|  Takes a set of data associated with a given cluster gets move probability
|  associated with the provided latent pattern
|
---------------------------------------------------------------------------------*/
MrBFlt ForceLatentPatternWithChangedData(int *oldAllocationVector, int *newAllocationVector, int *oldLatentMatrix, int *updatedLatentMatrix, int numChars, int compMatrixStart, int oldTable, int newTable, MrBFlt rho, RandLong *seed)
{
    int         i, numCharsInClusterSrc=0, latentIdxSrc=-1, *dataSrc, numMissingSrc,
                origLatentStatesSrc[numLocalTaxa], latentIdxDst=-1, *dataDst, numMissingDst,
                numCharsInClusterDst=0, origLatentStatesDst[numLocalTaxa];
    MrBFlt      backwardsMoveProbSrc, backwardsMoveProbDst, backwardsMoveProbTotal;

    // Get number of characters in source cluster and index of original latent pattern
    for (i=0; i<numChars; i++)
        if (oldAllocationVector[i] == oldTable)
        {
            numCharsInClusterSrc++;
            if (latentIdxSrc < 0)
                latentIdxSrc = i;
        }

    // Get data from characters that belong to the newTable cluster */
    dataSrc = GetClusterData(oldAllocationVector, oldTable, numCharsInClusterSrc, numChars, compMatrixStart);

    // Get number of missing characters
    numMissingSrc = GetNumMissingChars(numCharsInClusterSrc * numLocalTaxa, dataSrc);

    // Get latent pattern associated with this cluster
    for (i=0; i<numLocalTaxa; i++)
        origLatentStatesSrc[i] = oldLatentMatrix[pos(i, latentIdxSrc, numChars)];

    // Finally get the move prob associated with the original latent states
    backwardsMoveProbSrc = ForceLatentPattern(dataSrc, origLatentStatesSrc, numCharsInClusterSrc, numMissingSrc, rho, seed);

    //
    // printf("src data:\n");
    // for (i=0; i<numLocalTaxa; i++)
    // {
    //     for (int j=0; j<numCharsInClusterSrc; j++)
    //         printf("%c ", WhichStand(dataSrc[pos(i,j,numCharsInClusterSrc)]));
    //     printf("\n");
    // }
    // printf("\n");
    //
    // printf("orig latent states src: ");
    // for (i=0; i<numLocalTaxa; i++)
    //     printf("%c ",WhichLatent(origLatentStatesSrc[i]));
    // printf("\n");
    //
    // printf("backwards move prob src: %f\n",backwardsMoveProbSrc);

    // Free allocations
    free (dataSrc);
    dataSrc = NULL;


    // Get number of characters in destination cluster and index of destination latent pattern after move
    for (i=0; i<numChars; i++)
        if (oldAllocationVector[i] == newTable)
        {
            numCharsInClusterDst++;
            if (latentIdxDst < 0)
                latentIdxDst = i;
        }

    if (numCharsInClusterDst > 0) // It's possible that the cluster has disappeared
    {
        // Get data from characters that belong to the newTable cluster */
        dataDst = GetClusterData(oldAllocationVector, newTable, numCharsInClusterDst, numChars, compMatrixStart);

        // Get number of missing characters
        numMissingDst = GetNumMissingChars(numCharsInClusterDst * numLocalTaxa, dataDst);

        // Get latent pattern associated with this cluster
        for (i=0; i<numLocalTaxa; i++)
            origLatentStatesDst[i] = oldLatentMatrix[pos(i, latentIdxDst, numChars)];

        // Finally get the move prob associated with the original latent states
        backwardsMoveProbDst = ForceLatentPattern(dataDst, origLatentStatesDst, numCharsInClusterDst, numMissingDst, rho, seed);

        // printf("dst data:\n");
        // for (i=0; i<numLocalTaxa; i++)
        // {
        //     for (int j=0; j<numCharsInClusterDst; j++)
        //         printf("%c ", WhichStand(dataDst[pos(i,j,numCharsInClusterDst)]));
        //     printf("\n");
        // }
        // printf("\n");
        //
        // printf("orig latent states dst: ");
        // for (i=0; i<numLocalTaxa; i++)
        //     printf("%c ",WhichLatent(origLatentStatesDst[i]));
        // printf("\n");
        //
        // printf("backwards move prob dst: %f\n",backwardsMoveProbDst);

        // Free allocations
        free (dataDst);
        dataDst = NULL;
    }
    else
    {
        backwardsMoveProbDst = 1.0;
    }

    backwardsMoveProbTotal = log(backwardsMoveProbSrc) + log(backwardsMoveProbDst);

    return backwardsMoveProbTotal;
}



/*---------------------------------------------------------------------------------
|
|  GetNumMissingChars
|
|  Gets number of missing characters from a set of characters (with length = numChars)
|
---------------------------------------------------------------------------------*/
int GetNumMissingChars(int numChars, int *pattern)
{
    int         i, numMissing=0;

    for (i=0; i<=numChars; i++)
        {
        if ((pattern[i] == MISSING) || (pattern[i] == GAP))
            numMissing++;
        }

    return numMissing;
}


/*---------------------------------------------------------------------------------
|
|  PrintLatentMatrixToScreen
|
|  Prints latent matrix.
|
---------------------------------------------------------------------------------*/
int PrintLatentMatrixToScreen(int numChars, int *latentMatrix)
{
    int         i, j;

    MrBayesPrint("Latent matrix:\n");
    for (i=0; i<numLocalTaxa; i++)
        {
        MrBayesPrint("\t");
        for (j=0; j<numChars; j++)
            MrBayesPrint("%c ",WhichLatent(latentMatrix[pos(i,j,numChars)]));
        MrBayesPrint("\n");
        }
    MrBayesPrint("\n\n");

    MrBayesPrint ("Press return to continue\n");
    getchar();

    return NO_ERROR;
}
