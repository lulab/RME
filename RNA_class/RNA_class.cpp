// Example of RNA Class Usage.
//

//Provide example of free energy calculation.
//Program takes one parameters.  The first is the name of a ct file to read.
//The program then output is the folding free energy change to std out. 

#include "RNA.h"
#include "../src/structure.h" //get the definition of ct just to read the ct from disk
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[])
{

	RNA *rna;
	structure ct;
	int s;
	int i,j;
	char *string;
	int error;
	double free_energy;


	if (argc!=3) {
		std::cout << "Usage: RNA_class ct_file_input structure_#_for_energy_calculation\n";
		return 0;
	}

	//open a ct for structure input for this demonstration
	openct(&ct,argv[1]);

	//convert the input for structure number to integer
	s = atoi(argv[2]);

	//Now use the ct information to encode structural information in RNA

	//First, build a string to input sequence to RNA
	string = new char [ct.numofbases+1];
	for (i=0;i<=ct.numofbases;i++) string[i]=ct.nucs[i+1];
	string[ct.numofbases]='\0'; //string must be null terminated

	rna = new RNA(string); //construct instance of RNA

	//now encode structure information
	for (i=1;i<=ct.numofstructures;i++) {
		for (j=1;j<=ct.numofbases;j++) {

			//identify pairs
			if (ct.basepr[i][j]>j) {
			
				//This is how pairs are specified in RNA
				//j is paired to ct->basepr[i][j]
				//i is the structure number
				//i defaults to 1, so don't specify the structure number for a single structure case
				error = rna->SpecifyPair(j,ct.basepr[i][j],i);
				if(error!=0) {
					
					//check to make sure that the return is zero, or else an error has occured
					std::cerr << rna->GetErrorMessage(error);
					return 0;
					


				}

			}
		}
	}

	//now return the calulated energy and display it:
	free_energy = rna->CalculateFreeEnergy(s);
	error = rna->GetErrorCode();
	if (error==0) {
		//Note that when calculate energy is called the first time, RNA reads parameter files from
		//disk at the location specified by environment variable DATAPATH.

		//These are the .dat files found in the data_tables directory of RNAstructure
		std::cout << "Free energy change is: "<<rna->CalculateFreeEnergy(s) << "\n";
	}
	else {
		std::cerr << rna->GetErrorMessage(error);
	}
		


	delete rna;
	delete[] string;
	return 0;
}

