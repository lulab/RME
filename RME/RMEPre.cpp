#include <iostream>
#include <cstdlib>
using namespace std;
// from RNAstructure package
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"
#include "../RNA_class/RNA.h"
#include "../src/boltzmann.h"

#include "RMEPre.h"
#include "Utils.h"

RMEPre::RMEPre()
{
	InitParams();
}

RMEPre::RMEPre(const std::string& seqFile, const std::string& probFile)
{
	InitParams();
	Create(seqFile, probFile);
}

RMEPre::RMEPre(const std::string& seq, std::list<RMEProbRecord>& probRecords)
{
	InitParams();
	Create(seq, probRecords);
}

RMEPre::RMEPre(RNA* strand, std::list<RMEProbRecord>& probRecords)
{
	InitParams();
	Create(strand, probRecords);
}

RMEPre::~RMEPre()
{
	if(strand && !strandGiven)
		delete strand;
}

void RMEPre::Create(const std::string& seqFile, const std::string& probFile)
{
	strand = new RNA(seqFile.c_str(), 2, true);
	probRecords = ReadProbFile(probFile);
	created = true;
}

void RMEPre::Create(const std::string& seq, std::list<RMEProbRecord>& probRecords)
{
	strand = new RNA(seq.c_str(), true);
	this->probRecords = probRecords;
	created = true;
}

void RMEPre::Create(RNA* strand, std::list<RMEProbRecord>& probRecords)
{
	this->strand = strand;
	this->probRecords = probRecords;
	created = true;
	strandGiven = true;
}

void RMEPre::InitParams()
{
	m = 0.5;
	epsilon = 0.01;
	temperature = 310.25;
	created = false;
	strand = NULL;
	strandGiven = false;
}

bool RMEPre::CalcPartitionFunction(const std::string& pfsFile)
{
	if(!created)
		Die("The RMEPre objected is not created before calculating partition function");
	ErrorChecker<RNA>* checker = new ErrorChecker<RNA>(strand);
	int error = checker->isErrorStatus();
	if(error == 0 && (temperature != 310.15))
	{
		int tempError = strand->SetTemperature(temperature);
		error = checker->isErrorStatus(tempError);
	}
	if(error == 0)
	{
		//cout << "Begining calculating partition function." << endl;
		//TProgressDialog* progress = new TProgressDialog();
		//strand->SetProgress(*progress);
		CalculatePseudoEnergy(probRecords);
		int partError = strand->PartitionFunction(pfsFile.c_str());
		error = checker->isErrorStatus(partError);
		//cout << "Finishing calculating partition function." << endl;
		//strand->StopProgress();
		//delete progress;
	}
	delete checker;
	return (error == 0);
}

void RMEPre::CalculatePseudoEnergy(std::list<RMEProbRecord>& probRecords)
{
	structure* ct = strand->GetStructure();
	ct->shaped = true;
	ct->SHAPE = new double[2*ct->numofbases + 1];
	ct->SHAPEss = new double[2*ct->numofbases + 1];
	// default pseudo energy is 0
	for(int i = 1; i <= 2*ct->numofbases; i ++)
		ct->SHAPE[i] = 0;
	
	for(std::list<RMEProbRecord>::iterator iter = probRecords.begin(); iter != probRecords.end(); ++ iter)
	{
		double energy = CalculatePseudoEnergy(iter->value);
		ct->SHAPE[iter->pos] = energy;
		// duplicate the array
		ct->SHAPE[iter->pos + ct->numofbases] = energy;
	}
}

double RMEPre::CalculatePseudoEnergy(double data)
{
	double energy = 0;
	const double kT = 6.163341;
		
	if ((0 <= data) && (data <= 1)) {
		double w = (data + epsilon)/(1.0 - data + epsilon);
		energy = -kT*m*log(w);
	}
	else
		energy = 0;
	return energy;
}

void RMEPre::AddOptions(ParseCommandLine* parser)
{
	epsilonOptions.push_back( "-epsilon");
	epsilonOptions.push_back( "-Epsilon");
	epsilonOptions.push_back( "--Epsilon");
	parser->addOptionFlagsWithParameters( epsilonOptions, "A small value used when the constraint probability falls below 0. Default is 0.01.");

	mOptions.push_back( "-m");
	parser->addOptionFlagsWithParameters( mOptions, "Constraint weight. Default is 1.0.");
	
	tempOptions.push_back( "-t" );
	tempOptions.push_back( "-T" );
	tempOptions.push_back( "--temperature" );
	parser->addOptionFlagsWithParameters( tempOptions, "Specify the temperature at which calculation takes place in Kelvin. Default is 310.15 K, which is 37 degrees C." );
}

bool RMEPre::GetOptions(ParseCommandLine* parser)
{
	if(parser->contains(tempOptions))
	{
		parser->setOptionDouble(tempOptions, temperature);
		if(temperature < 0)
		{
			parser->setError("temperature");
			return false;
		}
	}
	if(parser->contains(epsilonOptions))
		parser->setOptionDouble(epsilonOptions, epsilon);
	if(parser->contains(mOptions))
		parser->setOptionDouble(mOptions, m);
	return true;
}
