#ifndef __RMEPOST_H__
#define __RMEPOST_H__

#include <string>
#include <vector>
#include <list>

#include "Utils.h"

// forward declarations
class RNA;
class structure;
class ParseCommandLine;

struct RNAHelix
{
	int start;
	int rstart;
	int end;
	int length;
	int score;
	int curScore;
	
	RNAHelix():
	start(-1), rstart(-1), length(0), score(0), curScore(0)
	{}
};

class RMEPost
{
public:
	// construct from files
	RMEPost();
	RMEPost(const std::string& pfsFile, 
				 const std::string& probFile);
	// construct from an RNA object with a given list of probability
	RMEPost(RNA* strand, std::list<RMEProbRecord>& probRecords);
	~RMEPost();
	// lazy creation
	void Create(const std::string& pfsFile, const std::string& probFile);
	void Create(RNA* strand, std::list<RMEProbRecord>& probRecords);
	
	void SetName(const std::string& name);
	// prediction
	void Fold();
	structure* GetStructure() { return this->ct; }
	// save structure to file
	void SaveAsDotBracket(const std::string& fileName);
	void SaveAsCt(const std::string& fileName);
	
	// for commandline parsing
	void AddOptions(ParseCommandLine* parser);
	// check and set parameters, returns true if successful
	bool GetOptions(ParseCommandLine* parser);
public:
	// parameters
	double gamma1;
	double gamma2;
	double gamma;
	int minHelixLen;
	int maxHelixLen;
private:
	void Init(RNA* strand);
	void InitParams();
	//update BPPM
	void UpdateBppm();

	//MaxExpect
	void MEPredict();
	void Fill(Matrix<double>& V, Matrix<double>& W, Matrix<int>& T);
	void TraceBack(const Matrix<int>& T);
private:
	bool created;
	structure* ct;
	std::string name;
	int length;	//length of the RNA sequence
	std::vector<int> enc;	//encoded RNA sequence
	//matrices and vectors
	Matrix<double> bppm;	//base pair probability matrix
	std::vector<double> bpp;	//base pair probability vector
	std::vector<double> q;	//constraint vector
	std::list<RMEProbRecord> probRecords;	//records in the probability file
	Matrix<short> delta;	//indicator matrix, 0, 1 for unpaired and paired
	
	// for command line parsing
	std::vector<std::string> gamma1Options;
	std::vector<std::string> gamma2Options;
	std::vector<std::string> gammaOptions;
	std::vector<std::string> minHelixLenOptions;
	std::vector<std::string> maxHelixLenOptions;
};

#endif