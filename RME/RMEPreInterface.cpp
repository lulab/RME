#include <iostream>
#include <cstdlib>
using namespace std;
// from RNAstructure package
#include "../src/ParseCommandLine.h"
#include "../RNA_class/RNA.h"

#include "RMEPre.h"
#include "Utils.h"

int main(int argc, char** argv)
{
	RMEPre pre;
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine("RME-Pre");
	// required arguments
	parser->addParameterDescription("seq file", "The name of a file containing an input sequence.");
	parser->addParameterDescription("prob file", "The name of a file containing probability from structural probing data." );
	parser->addParameterDescription("pfs file", "The name of a partition function save file to be write to." );
	// add RME-Post options
	pre.AddOptions(parser);
	// parser command line
	parser->parseLine(argc, argv);
	// get required parameters from command line arguments
	string seqFile = parser->getParameter(1);
	string probFile = parser->getParameter(2);
	string pfsFile = parser->getParameter(3);
	// get RME-Pre options
	if(pre.GetOptions(parser))
	{
		pre.Create(seqFile, probFile);
		pre.CalcPartitionFunction(pfsFile);
	}
	else
		return -1;
	return 0;
}