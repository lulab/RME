/*
 * A program that converts a dot bracket file to a CT file.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "dot2ct.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
dot2ct_Interface::dot2ct_Interface() {

	// Initialize the calculation type description.
	calcType = "Dot bracket file conversion";
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool dot2ct_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "dot2ct" );
	parser->addParameterDescription( "bracket file", "The name of a file containing the dot bracket structure to convert." );
	parser->addParameterDescription( "ct file", "The name of a CT file to which output will be written." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		bracketFile = parser->getParameter( 1 );
		ctFile = parser->getParameter( 2 );
	}

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void dot2ct_Interface::run() {

	// Show a message saying that conversion has started.
	cout << "Converting dot bracket file..." << flush;

	// Create a variable that handles errors.
	int error = 0;

	// Define the maximum length of the line in the dot bracket file.
	const int MAXLINELENGTH = 15000;

	// Define variables needed for converstion.
	ifstream in;
	char line[MAXLINELENGTH];
	structure ct;
	bool *right;
	int i,j;

	// Open the .bracket file.
	in.open( bracketFile.c_str() );

	// Copy the first line to the CT label and make sure it has a newline on the end of it.
	in.getline( line,MAXLINELENGTH );
	strcpy( ct.ctlabel[1], line + 1 );
	strcat( ct.ctlabel[1], "\n" );

	// Get the second line and use it to allocate the appropriate amount of space for the structure.
	// There's always one structure only, at the current time.
	in.getline( line,MAXLINELENGTH );
	ct.numofbases = strlen( line );
	ct.allocate( ct.numofbases + 1 );
	ct.numofstructures = 1;

	// Read in the nucleotides.
	for( i = 0; i < strlen( line ); i++ )
		ct.nucs[i+1] = line[i];

	in.getline( line,MAXLINELENGTH );

	// Initialize the pairs to zero.
	for( i = 1; i <= ct.numofbases; i++ ) ct.basepr[1][i] = 0;

	// Parse the dots and brackets.
	right = new bool[ strlen( line ) ];
	for( i=0; i < strlen( line ); i++ ) {
		if( line[i] == '(' || line[i] == '>' ) right[i] = true;
		else right[i] = false;
	}

	// Find left-facing brackets.
	for( i = 0; i < strlen( line ); i++ ) {
		if( line[i] == ')' || line[i] == '<' ) {

			// Find a 3' paired nuc, then find the 5' partner.
			j = i - 1;
			while( (right[j] == false) && j >= 0 ) { j--;}

			right[j]=false;

			// Register the pair that was found.
			ct.basepr[1][i+1] = j + 1;
			ct.basepr[1][j+1] = i + 1;
		}
	}

	// Before outputting the ct file, the historical numbering needs to be set.
	for( i = 1; i <= ct.numofbases; i++ ) ct.hnumber[i] = i;
	ctout( &ct, ctFile.c_str() );

	// Close the file and delete the bracket data array.
	in.close();
	delete[] right;

	// Show a message saying conversion is done.
	cout << "done." << endl;

	// Print confirmation of run finishing.
	if( error == 0 ) { cout << calcType << " complete." << endl; }
	else { cerr << calcType << " complete with errors." << endl; }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	dot2ct_Interface* runner = new dot2ct_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
