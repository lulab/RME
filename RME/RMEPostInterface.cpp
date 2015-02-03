#include <iostream>
#include <cstdlib>
using namespace std;
// from RNAstructure package
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"
#include "../RNA_class/RNA.h"

#include "RMEPost.h"
#include "Utils.h"

int main(int argc, char** argv)
{
	RMEPost post;
	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine("RME-Post");
	// required arguments
	parser->addParameterDescription("pfs file", "The name of a partition function save file from partition-cons." );
	parser->addParameterDescription("prob file", "Probability from structural probing data." );
	parser->addParameterDescription("ct file", "The output structure file in ct format." );
	// add RME-Post options
	post.AddOptions(parser);
	// parser command line
	parser->parseLine(argc, argv);
	
	// get required parameters from command line arguments
	string pfsFile = parser->getParameter(1);
	string probFile = parser->getParameter(2);
	string ctFile = parser->getParameter(3);
	// get RME-Post options
	if(post.GetOptions(parser))
	{
		post.Create(pfsFile, probFile);
		post.Fold();
		post.SaveAsCt(ctFile);
	}
	else
		return -1;
	return 0;
}