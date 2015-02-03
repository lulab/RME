/*
 * A program that converts a CT file to a dot bracket file.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "ct2dot.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ct2dot_Interface::ct2dot_Interface() {

	// Initialize the calculation type description.
	calcType = "CT file conversion";
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool ct2dot_Interface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ct2dot" );
	parser->addParameterDescription( "ct file", "The name of a file containing the CT structure to convert." );
	parser->addParameterDescription( "structure number", "The number, one-indexed, of the structure to convert in the CT file." );
	parser->addParameterDescription( "bracket file", "The name of a dot bracket file to which output will be written." );

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	string numberString;
	if( !parser->isError() ) {
		ctFile = parser->getParameter( 1 );
		numberString = parser->getParameter( 2 );
		bracketFile = parser->getParameter( 3 );
	}

	// Convert the structure number parameter into an integer.
	// If that can't be done, set an error.
	int testInt = 0;
	stringstream testStream( numberString );
	if( testStream >> testInt ) { number = testInt; }
	else { parser->setError( "structure number" ); }

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void ct2dot_Interface::run() {

	// Show a message saying that conversion has started.
	cout << "Converting CT file..." << flush;

	// Create a variable that handles errors.
	int error = 0;

	// Initialize and open the CT file to convert.
	structure ct;
	openct( &ct, ctFile.c_str() );

	// Initialize and open the output dot bracket file.
	ofstream out;
	out.open( bracketFile.c_str() );

	// Write the comment line.
	// Make sure that the line ends in a newline, if not, add a newline.
	out << ">" << ct.ctlabel[number];
	if( ct.ctlabel[number][strlen(ct.ctlabel[number])-1] != '\n' ) {
		out << "\n";
	}

	// Write the sequence.
	for( int i = 1; i <= ct.numofbases; i++ ) {
		out << ct.nucs[i];
	}
	out << "\n";

	// Write the conversion mask.
	for( int i = 1; i <= ct.numofbases; i++ ) {
		if( ct.basepr[number][i] > i ) { out << "("; }
		else if( ct.basepr[number][i] == 0 ) { out << "."; }
		else { out << ")"; }
	}
	out << "\n";

	// Close the dot bracket file.
	out.close();

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

	ct2dot_Interface* runner = new ct2dot_Interface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
