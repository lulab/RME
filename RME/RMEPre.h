#ifndef __RMEPRE_H__
#define __RMEPRE_H__

#include <string>
#include <list>

// forward declarations
class RNA;
class ParseCommandLine;
struct RMEProbRecord;

class RMEPre
{
public:
	RMEPre();
	RMEPre(const std::string& seqFile, const std::string& probFile);
	RMEPre(const std::string& seq, std::list<RMEProbRecord>& probRecords);
	RMEPre(RNA* strand, std::list<RMEProbRecord>& probRecords);
	~RMEPre();
	// lazy creation
	void Create(const std::string& seqFile, const std::string& probFile);
	void Create(const std::string& seq, std::list<RMEProbRecord>& probRecords);
	void Create(RNA* strand, std::list<RMEProbRecord>& probRecords);
	
	bool CalcPartitionFunction(const std::string& pfsFile = "");
	void SaveAs(const std::string& pfsFile);
	RNA* GetStrand() { return strand; }
	
	// for commandline parsing
	void AddOptions(ParseCommandLine* parser);
	// check and set parameters, returns true if successful
	bool GetOptions(ParseCommandLine* parser);
public:
	// parameters
	double m;
	double epsilon;
	double temperature;
private:
	void InitParams();
	// calculate pseudo energy for all bases
	void CalculatePseudoEnergy(std::list<RMEProbRecord>& probRecords);
	// calculate pseudo energy from probability
	double CalculatePseudoEnergy(double data);
private:
	RNA* strand;
	bool strandGiven;
	bool created;
	std::list<RMEProbRecord> probRecords;
	// for commandline parsing
	std::vector<std::string> epsilonOptions;
	std::vector<std::string> mOptions;
	std::vector<std::string> tempOptions;
};

#endif