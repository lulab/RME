#ifndef __RME_H__
#define __RME_H__

#include <string>
#include <vector>

#define MAX_PATH_LENGTH		4096
#define MAX_SEQ_LENGTH		5000
#define MIN_LOOP_LENGTH		4

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
	RMEPost(const std::string& pfsFile, 
				 const std::string& probFile);
	// construct from an RNA object with partition function calculated
	RMEPost(RNA* strand, const std::string& probFile);
	// construct from an RNA object with a given list of probability
	RMEPost(RNA* strand, std::list<RMEProbRecord>& probRecords);
	~RMEPost();
	
	void Fold();
	void SaveAsDotBracket(const std::string& fileName);
	void SaveAsCt(const std::string& fileName);
	
	// Parameter settings
	void SetGamma1(double value) {this->gamma1 = value;}
	void SetGamma2(double value) {this->gamma2 = value;}
	void SetGamma(double value) {this->gamma = value;}
	void SetMinHelixLen(int value) {this->minHelixLen = value;}
	void SetMaxHelixLen(int value) {this->maxHelixLen = value;}
private:
	void Init(RNA* strand);
	//update BPPM
	void UpdateBppm();

	//MaxExpect
	void MEPredict();
	void Fill(Matrix<double>& V, Matrix<double>& W, Matrix<int>& T)
	void TraceBack(const Matrix<int>& T)
private:
	// parameters
	double gamma1;
	double gamma2;
	double gamma;
	int minHelixLen;
	int maxHelixLen;
	
	structure* ct;
	bool strandGiven; 
	int length;	//length of the RNA sequence
	std::vector<int> enc;	//encoded RNA sequence
	//matrices and vectors
	Matrix<double> bppm;	//base pair probability matrix
	std::vector<double> bpp;	//base pair probability vector
	std::vector<double> ssp;	//single-stranded probability vector
	std::vector<double> q;	//constraint vector
	std::list<RMEProbRecord> probRecords;	//records in the probability file
	Matrix<short> delta;	//indicator matrix, 0, 1 for unpaired and paired
};

struct RMEPreInterface
{
	// parameters
	double m;
	double epsilon;;
	double temperature;
	// parser options
	std::vector<std::string> epsilonOptions;
	std::vector<std::string> mOptions;
	std::vector<std::string> tempOptions;
	
	RMEPreInterface();
	void AddOptions(ParseCommandLine* parser);
	int GetOptions(ParseCommandLine* parser);
	bool ReadProb(const string& probFile);
	// calculate pseudo energy for all bases
	void CalculatePseudoEnergy(std::list<RMEProbRecord>& probRecords);
	// calculate pseudo energy from probability
	void CalculatePseudoEnergy(double data);
	void Run(RNA* strand, std::list<RMEProbRecord>& probRecords);
	void Run(const std::string& seqFile, const std::string& probFile, const std::string& pfsFile);
	void Run(RNA* strand, const std::string& probFile, const std::string& pfsFile);
	
};

struct RMEPostInterface
{
	// parameters
	double gamma1;
	double gamma2;
	double gamma;
	int minHelixLen;
	int maxHelixLen;
	// parser options
	std::vector<std::string> gamma1Options;
	std::vector<std::string> gamma2Options;
	std::vector<std::string> gammaOptions;
	std::vector<std::string> minHelixLenOptions;
	std::vector<std::string> maxHelixLenOptions;
	
	RMEPostInterface();
	void AddOptions(ParseCommandLine* parser);
	int GetOptions(ParseCommandLine* parser);
	void Run(const std::string& pfsFile, const std::string& probFile, const std::string& ctFile);
	void Run(RNA* strand, const std::string& probFile, const std::string& ctFile);
};



#endif