//Calculate scores when comparing a predicted structure (called test) to a known structure (called correct).
#include "score.h"


//Calculate PPV, the fraction of predicted pairs that are correct.
	//Requires pointers to two structure classes, the test and correct structures
	//score is a pointer to integer that tabulates the number of correct predictions
	//basepairs is the total number of predicted pairs
	//structure number is the number of the structure to score in test (correct can only have one structure)
	//exact is a bool that indicates whether exact matches are required, normally this is set false to allow slippage

//This function returns an int that prevides an error code
	//0 = no error
	//1 = more than one structure in correct
	//2 = two structures are different lengths
	
int scorerppv(const structure* correct, const structure* test, int *score, int *basepairs, const int structurenumber, const bool exact){
	int i,j;

	//initialize the pair counts
	(*basepairs)=0;
	*score=0;

	//check for errors
	if (correct->numofstructures!=1) {
		
		return 1;
	}
	if (correct->numofbases!=test->numofbases) {
		
		return 2;
	}



	//count the bases in the predicted! ct:
	for (i=1;i<=correct->numofbases;i++) {
		if (test->basepr[structurenumber][i]>i) (*basepairs)++;
	}


	j = structurenumber;

	
	
	for (i=1;i<=test->numofbases;i++) {
		
		if (test->basepr[j][i]>i) {
      			if (test->basepr[j][i]==correct->basepr[1][i])
         			(*score)++;
         		else if (!exact&&(test->basepr[j][i])==correct->basepr[1][i]+1)
         			(*score)++;
         		else if (!exact&&(test->basepr[j][i])==correct->basepr[1][i]-1)
         			(*score)++;
				else if (!exact&&(test->basepr[j][i])==correct->basepr[1][i+1])
					(*score)++;
				else if (!exact&&(test->basepr[j][i])==correct->basepr[1][i-1])
					(*score)++;
				
		}
   	}

	
	return 0;

}

//Calculate sensitivity, the fraction of known pairs correctly predicted.
	//Requires pointers to two structure classes, the test and correct structures
	//score is a pointer to integer that tabulates the number of correct predictions
	//basepairs is the total number of known pairs
	//structure number is the number of the structure to score in test (correct can only have one structure)
	//exact is a bool that indicates whether exact matches are required, normally this is set false to allow slippage

//This function returns an int that prevides an error code
	//0 = no error
	//1 = more than one structure in correct
	//2 = two structures are different lengths

int scorer(const structure* correct, const structure* test, int *score, int *basepairs, const int structurenumber, const bool exact){
	int i,j;

	//initialize the pair counts
	(*basepairs)=0;
	*score=0;

	if (correct->numofstructures!=1) {
		return 1;
	}
	if (correct->numofbases!=test->numofbases) {
		return 2;
	}


	
	//count the bases in the phylogenetic ct:
	for (i=1;i<=correct->numofbases;i++) {
		if (correct->basepr[1][i]>i) (*basepairs)++;
	}

	//Now check every pair

	j = structurenumber;	
	
	for (i=1;i<=test->numofbases;i++) {
		
    		if (correct->basepr[1][i]>i) {
      			if (test->basepr[j][i]==correct->basepr[1][i])
         			(*score)++;
         		else if (!exact&&(test->basepr[j][i]+1)==correct->basepr[1][i])
         			(*score)++;
         		else if (!exact&&(test->basepr[j][i]-1)==correct->basepr[1][i])
         			(*score)++;
				else if (!exact&&(test->basepr[j][i+1])==correct->basepr[1][i])
					(*score)++;
				else if (!exact&&(test->basepr[j][i-1])==correct->basepr[1][i])
					(*score)++;
				
      		}
   	}

	return 0;

}
