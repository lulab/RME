#include <cmath>
#include <cstdlib>
#include <stdarg.h>
#include <cctype>
#include <iostream>
#include <sstream>
#include <map>
using namespace std;
// system calls
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>

// from RNAstructure package
#include <src/score.h>

#include "Utils.h"

void Die(const char* format, ...)
{
	char buffer[256];
	va_list args;
	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
	cerr << "Fatal error: " << buffer << endl;
	exit(-1);
}

std::vector<int> EncodeRNA(const std::string& seq)
{
	std::vector<int> enc(seq.size());
	for(size_t i = 0; i < seq.size(); i ++)
	{
		bool valid = false;
		char base = toupper(seq[i]);
		if(base == 'T')
			base = 'U';
		for(int j = 0; j < 4; j ++)
		{
			if(base == BASE[j])
			{
				enc[i] = j;
				valid = true;
				break;
			}
		}
		if(!valid)
			Die("invalid character '%c' found in RNA seqeunce", seq[i]);
	}
	return enc;
}


//Check if a file exists
bool FileExists(const std::string& fileName)
{
	ifstream fin(fileName.c_str());
	if(fin)
		return true;
	else
		return false;
}

std::list<RMEProbRecord> ReadProbFile(const std::string& fileName)
{
	ifstream fin(fileName.c_str());
	if(!fin)
		Die("cannot read the probability file: ", fileName.c_str());
	
	std::list<RMEProbRecord> records;
	
	while(!fin.eof())
	{
		RMEProbRecord record;
		fin >> record.pos >> record.value;
		records.push_back(record);
	}
	fin.close();
	return records;
}

std::vector<double> ReadProbFile(const std::string& fileName, int length)
{
	vector<double> prob(length, -1);
	
	ifstream fin(fileName.c_str());
	if(!fin)
		Die("cannot read the probability file: ", fileName.c_str());
	
	int lineno = 1;
	while(!fin.eof())
	{
		string line;
		std::getline(fin, line);
		int pos;
		double value;
		
		istringstream is(line);
		is >> pos >> value;
		if((pos < 1) || (pos > length))
			Die("Index not between 1 and %d in the probabiity file at line %d", length, lineno);
		
		prob[pos - 1] = value;
		lineno ++;
	}
	fin.close();
	
	return prob;
}

std::string GetDataPath(const char* argv)
{
	string progName(argv);
	size_t i = progName.rfind('/');
	if(i == std::string::npos)
		return string("/../data_tables/");
	else
		return progName.substr(0, i) + "/../data_tables/";
}

// scores are written to sensitivity and ppv
void GetScore(structure* predct, structure* refct, double& sensitivity, double& ppv)
{
	int positive = 0;
	int pairs = 0;
	bool exact = false;
	
	scorer(refct, predct, &positive, &pairs, 1, exact);
	sensitivity = (double(positive)/double(pairs))*100;
	scorerppv(refct, predct, &positive, &pairs, 1, exact);
	if(pairs <= 0)
		ppv = 0.0;
	else
		ppv = (double(positive)/double(pairs))*100;
}

std::vector<std::string> SplitString(const std::string& s, char delim)
{
	size_t start = 0;
	std::vector<std::string> tokens;
	for(size_t i = 0; i < s.size(); i ++)
	{
		if(s[i] == delim)
		{
			tokens.push_back(s.substr(start, i - start));
			start = i + 1;
		}
	}
	if(start != s.size())
		tokens.push_back(s.substr(start, s.size() - start));
	return tokens;
}

bool MakeDir(const std::string& fileName)
{
#ifdef WIN32
	int status = mkdir(fileName.c_str());
#else
	int status = mkdir(fileName.c_str(), 0777);
#endif
	return (status == 0);
}

bool DirExists(const std::string& dirName)
{
	struct stat statbuf;
	if(stat(dirName.c_str(), &statbuf) == 0)
		return true;

	//cerr << "Cannot stat " << dirName << ": " << strerror(errno) << endl;
	return false;
}

std::vector<RMEInputDataRecord> ReadRMEInputData(const std::string& fileName, bool hasRef)
{
	// read training data
	ifstream fin(fileName.c_str());
	if(!fin)
		Die("Cannot read input data: %s", fileName.c_str());
	map<string, RMEInputDataRecord> dataMap;
	// skip header
	string line;
	std::getline(fin, line);
	int lineno = 2;
	while(!fin.eof())
	{
		std::getline(fin, line);
		istringstream is(line);
		string name;
		int index;
		string prob;
		char base;
		int refbp;

		if(hasRef)
			is >> name >> index >> prob >> base >> refbp;
		else
			is >> name >> index >> prob >> base;
		
		if(is.bad())
			Die("Failed to read input data at line %d", lineno);
		if(name == "")
			break;
		if(dataMap.find(name) == dataMap.end())
		{
			// create a new record
			RMEInputDataRecord newRecord;
			newRecord.seq.resize(MAX_SEQ_LENGTH, 'N');
			newRecord.length = 0;
			newRecord.name = name;
			if(hasRef)
				newRecord.refbp.resize(MAX_SEQ_LENGTH, -1);
			
			dataMap[name] = newRecord;
		}
		RMEInputDataRecord& record = dataMap[name];
		if(index < 1 || index > MAX_SEQ_LENGTH)
			Die("The index in input data is not between 1 and %d at line %d", MAX_SEQ_LENGTH, lineno);
		if(index > record.length)
			record.length = index;
		if(prob != "NA")
		{
			RMEProbRecord probRecord;
			probRecord.pos = index;
			probRecord.value = strtod(prob.c_str(), NULL);
			record.probRecords.push_back(probRecord);
		}
		record.seq[index - 1] = base;
		if(hasRef)
			record.refbp[index - 1] = refbp;

		lineno ++;
	}
	fin.close();
	// check for missing data
	bool hasMissing = false;
	for(std::map<string, RMEInputDataRecord>::iterator iter = dataMap.begin();
	        iter != dataMap.end(); ++ iter)
	{
		string name = iter->first;
		RMEInputDataRecord& record = iter->second;
		record.seq.resize(record.length);
		if(hasRef)
			record.refbp.resize(record.length);
		for(int i = 0; i < record.length; i ++)
		{
			if(record.seq[i] == 'N')
			{
				hasMissing = true;
				cerr << name << " has missing data at position " << i + 1 << endl;
				break;
			}
		}
	}
	if(hasMissing)
		Die("Failed due to missing data in the input data");

	// convert data map to a vector
	vector<RMEInputDataRecord> inputData;
	for(std::map<string, RMEInputDataRecord>::iterator iter = dataMap.begin();
	        iter != dataMap.end(); ++ iter)
		inputData.push_back(iter->second);
	
	return inputData;
}