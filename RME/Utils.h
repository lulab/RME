#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <pthread.h>
// from RNAstructure package
#include <src/structure.h>


#define MAX_PATH_LENGTH		4096
#define MAX_SEQ_LENGTH		5000
#define MIN_LOOP_LENGTH		4

//Utilities for operation on RNA sequence
const char BASE[4] = {'A', 'U', 'C', 'G'};

const int BASEPAIRS[4][4] = {
	{0, 1, 0, 0},
	{1, 0, 0, 1},
	{0, 0, 0, 1},
	{0, 1, 1, 0}
};

//Check if a and b can form a cannonical pair
inline bool IsBasePair(int a, int b)
{
	return BASEPAIRS[a][b];
}

//Convert characters to integers
std::vector<int> EncodeRNA(const std::string& seq);

//Print error message and terminate the program
void Die(const char* format, ...);

//Check if a file exists
bool FileExists(const std::string& fileName);

bool DirExists(const std::string& dirName);
bool MakeDir(const std::string& fileName);

template <class T>
class Matrix
{
public:
	Matrix()
	{
		this->size1 = 0;
		this->size2 = 0;
	}
	
	Matrix(size_t size1, size_t size2)
	{
		resize(size1, size2);
	};
	
	Matrix(size_t size)
	{
		resize(size);
	}
	
	T operator()(size_t i, size_t j) const
	{
		return _data[i*size2 + j];
	}
	
	T& operator()(size_t i, size_t j)
	{
		return _data[i*size2 + j];
	}
	
	void resize(size_t size1, size_t size2)
	{
		this->size1 = size1;
		this->size2 = size2;
		_data.resize(size1*size2);
	};
	
	void resize(size_t size)
	{
		this->size1 = size;
		this->size2 = size;
		_data.resize(size1*size2);
	}
private:
	size_t size1;
	size_t size2;
	std::vector<T> _data;
};
// a record in the probability file
struct RMEProbRecord
{
	int pos;	// 1-based position
	double value;	// values between 0 and 1
};

// reads a probability file and returns a returns a list of records
std::list<RMEProbRecord> ReadProbFile(const std::string& fileName);
std::vector<double> ReadProbFile(const std::string& fileName, int length);

std::string GetDataPath(const char* argv);

// scores are written to sensitivity and ppv
void GetScore(structure* predct, structure* refct, double& sensitivity, double& ppv);

// split a string into a vector of tokens by a delimiter
std::vector<std::string> SplitString(const std::string& s, char delim);

// progress bar
template <class ValType>
class TProgressBar
{
public:
	TProgressBar(ValType maxVal, ValType step): curVal(ValType(0)), maxVal(maxVal), step(step) {}

	void Show()
	{
		cout << "\r";
		cout.flush();
		if(caption.size() > 1)
			cout << caption << ": ";
		std::streamsize prec = cout.precision();
		cout.precision(2);
		cout << std::fixed << double(curVal)/double(maxVal)*100.0
		     << " % (" << curVal << "/" << maxVal << ")";
		cout.precision(prec);
		cout.flush();
	}

	void SetCaption(const std::string caption) {
		this->caption = caption;
	}
	void Increment() {
		curVal += step;
	}
	void Set(ValType val) {
		curVal = (val < maxVal)? val : maxVal;
	}
	void Reset() {
		curVal = ValType(0);
	}
	void Reset(ValType maxVal, ValType step)
	{
		this->maxVal = maxVal;
		this->step = step;
		this->curVal = ValType(0);
	}
private:
	std::string caption;
	ValType curVal;
	ValType maxVal;
	ValType step;
};

// RME input data record
struct RMEInputDataRecord
{
	string name;
	string seq;
	list<RMEProbRecord> probRecords;
	int length;
	// reference structure is optional
	vector<int> refbp;
};

std::vector<RMEInputDataRecord> ReadRMEInputData(const std::string& fileName, bool hasRef = false);
// multi-threading functions
template <class ArgType>
struct ApplyData
{
	ArgType* arg;
	void (*func)(void*);
};


template <class ArgType>
void* ApplyThread(void* arg)
{
	static pthread_mutex_t mutex_q = PTHREAD_MUTEX_INITIALIZER;
	int retval = 0;
	while(true)
	{
		std::list<ApplyData<ArgType> >* q = (std::list<ApplyData<ArgType> >*)arg;
		pthread_mutex_lock(&mutex_q);
		if(q->size() < 1)
		{
			pthread_mutex_unlock(&mutex_q);
			pthread_exit(&retval);
		}
		ApplyData<ArgType> data = q->front();
		q->pop_front();
		pthread_mutex_unlock(&mutex_q);
		
		data.func(data.arg);
	}
	return NULL;
}


template <class ArgType>
void ApplyFunction(void (*func)(void*), std::vector<ArgType>& args, int nThreads)
{
	std::list<ApplyData<ArgType> > data;
	for(size_t i = 0; i < args.size(); i ++)
	{
		ApplyData<ArgType> item;
		item.arg = &(args[i]);
		item.func = func;
		data.push_back(item);
	}

	std::vector<pthread_t> threads(nThreads);
	for(int i = 0; i < nThreads; i ++)
		pthread_create(&(threads[i]), NULL, ApplyThread<ArgType>, (void*)&data);
	for(int i = 0; i < nThreads; i ++)
		pthread_join(threads[i], NULL);
}

#endif