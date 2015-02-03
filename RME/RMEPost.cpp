#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <stack>
using namespace std;
//from RNAstructure
#include "../RNA_class/RNA.h"
#include "../src/ParseCommandLine.h"

#include "Utils.h"
#include "RMEPost.h"

//// Implementation of class RMEPost ////
RMEPost::RMEPost()
{
	InitParams();
}

RMEPost::RMEPost(const std::string& pfsFile, 
				 const std::string& probFile)
{
	InitParams();
	Create(pfsFile, probFile);
}

RMEPost::RMEPost(RNA* strand, std::list<RMEProbRecord>& probRecords)
{
	InitParams();
	Create(strand, probRecords);
}

RMEPost::~RMEPost()
{
	if(ct)
		delete ct;
}

void RMEPost::Create(const std::string& pfsFile, const std::string& probFile)
{
	// read pfs file
	if(!FileExists(pfsFile.c_str()))
		Die("Cannot open the file: %s", pfsFile.c_str());
	RNA* strand = new RNA(pfsFile.c_str(), 3, true);
	probRecords = ReadProbFile(probFile);
	Init(strand);
	delete strand;
	created = true;
}

void RMEPost::Create(RNA* strand, std::list<RMEProbRecord>& probRecords)
{
	this->probRecords = probRecords;
	Init(strand);
	created = true;
}

void RMEPost::SetName(const std::string& name)
{
	this->name = name;
}

void RMEPost::InitParams()
{
	gamma1 = 0.0;
	gamma2 = 0.0;
	gamma = 0.5;
	minHelixLen = 0;
	maxHelixLen = 4;
	ct = NULL;
	length = 0;
	created = false;
}

void RMEPost::Init(RNA* strand)
{
	length = strand->GetSequenceLength();
	if(length <= 0)
		Die("RNA length <= 0");

	structure* in_ct = strand->GetStructure();
	//encode the RNA sequence
	string seq(length, 0);
	std::copy(in_ct->nucs + 1, in_ct->nucs + length + 1, seq.begin());
	enc = EncodeRNA(seq);
	//cout << seq << endl;
	//for(size_t i = 0; i < enc.size(); i ++)
	//	cout << enc[i];
	//cout << endl;
	// read BPPM
	//if(!in_ct->partitionfunctionallocated)
	//	Die("No partition function data available")
	bppm.resize(length);
	for(int i = 0; i < length; i ++)
	{
		for(int j = 0; j < length; j ++)
		{
			if(i < j)
				bppm(i, j) = strand->GetPairProbability(i + 1, j + 1);
			else
				bppm(i, j) = strand->GetPairProbability(j + 1, i + 1);
		}
	}
	// calculate marginal probability vector
	bpp.resize(length);
	for(int i = 0; i < length; i ++)
	{
		bpp[i] = 0.0;
		for(int j = 0; j < length; j ++)
			bpp[i] += bppm(i, j);
	}
	// initialize to marginal probability
	q = bpp;
	// read probability
	for(std::list<RMEProbRecord>::iterator iter = probRecords.begin(); iter != probRecords.end(); ++ iter)
		q[iter->pos - 1] = iter->value;
	
	//calculate delta matrix, which is a base pairing table
	delta.resize(length);
	for(int i = 0; i < length; i ++)
	{
		for(int j = i; j < length; j ++)
		{
			if(((j - i + 1) >= MIN_LOOP_LENGTH + 2) 
				&& (IsBasePair(enc[i], enc[j])))
				delta(i, j) = true;
			else
				delta(i, j) = false;
			delta(j, i) = delta(i, j);
		}
	}
	// make a copy of the structure
	ct = new structure(2);
	ct->numofbases = in_ct->numofbases;
	ct->allocate(ct->numofbases);
	for(int i = 1; i <= ct->numofbases; i ++)
	{
		ct->nucs[i] = in_ct->nucs[i];
		ct->hnumber[i] = in_ct->hnumber[i];
	}
	// initalize structural information
	for(int i = 1; i < ct->allocatedstructures; i ++)
		for(int j = 1; j <= ct->numofbases; j ++)
			ct->basepr[i][j] = 0;
	ct->numofstructures = 0;
}

void RMEPost::Fold()
{
	if(!created)
		Die("The RMEPost object has not been created before prediction");
	/*
	cout << "[RME-Post] gamma1=" << gamma1 << ", gamma2=" << gamma2 << endl;
	for(std::list<RMEProbRecord>::iterator iter = probRecords.begin(); 
		iter != probRecords.end(); ++ iter)
		cout << "q[" << iter->pos << "] = " << iter->value << endl;
	cout << "BPP before updating" << endl;
	for(size_t i = 0; i < bpp.size(); i ++)
		cout << bpp[i] << " ";
	cout << endl;
	*/
	UpdateBppm();
	/*
	cout << "BPP after updating" << endl;
	for(size_t i = 0; i < bpp.size(); i ++)
		cout << bpp[i] << " ";
	cout << endl;
	*/
	MEPredict();
}

void RMEPost::UpdateBppm()
{
	for(int t = 0; t < 2*length; t ++)
	{
		int i = (t > length)? (t - length) : 0;
		int j = (t < length)? t : length - 1;
		RNAHelix helix;
		while(i <= j)
		{
			if(!delta(i, j) || (i >= j))
			{
				//end of a helix
				if(helix.length > 0)
				{
					double w;
					if(helix.length < minHelixLen)
						w = 0.0;
					else if(helix.length >= maxHelixLen)
						w = 1.0;
					else
						w = double(helix.length)/(maxHelixLen);
					for(int l = 0; l < helix.length; l ++)
					{
						int i1 = helix.start + l;
						int j1 = helix.rstart - l;
						w *= gamma1;

						bppm(i1, j1) = (1.0 - w)*bppm(i1, j1) + w*q[i1]*q[j1];
						bppm(j1, i1) = bppm(i1, j1);
					}
				}
				helix.length = 0;
			}
			else
			{
				//start of a helix
				if(helix.length <= 0)
				{
					helix.start = i;
					helix.rstart = j;
				}
				//extend a helix
				helix.length ++;
			}
			i ++;
			j --;
		}
	}


	
	//update BPP
	for(int i = 0; i < length; i ++)
		bpp[i] = (1.0 - gamma2)*bpp[i] + gamma2*q[i];
}

/// MaxExpect algorithm
// traceback flags
#define TR_RIGHT		-1
#define TR_DOWN			-2
#define TR_PAIR			-3
#define TR_SS			-10

struct BasePair
{
	int i;
	int j;
	
	BasePair(int ii = -1, int jj = -1)
	: i(ii), j(jj) {}
};

void RMEPost::MEPredict()
{
	Matrix<double> V(length);
	Matrix<double> W(length);
	Matrix<int> T(length);
	
	ct->numofstructures = 1;
	ct->checknumberofstructures();
	ct->energy[1] = 0;
	strcpy(ct->ctlabel[1], name.c_str());
	
	Fill(V, W, T);
	TraceBack(T);
}

void RMEPost::Fill(Matrix<double>& V, Matrix<double>& W, Matrix<int>& T)
{
	//calculate single-stranded probability
	vector<double> ssp(length);
	for(int i = 0; i < length; i ++)
		ssp[i] = 1.0 - bpp[i];
	for(int l = 1; l <= length; l ++)
	{
		if(l <= MIN_LOOP_LENGTH)
		{
			for(int i = 0; i <= length - l; i ++)
			{
				int j = l - 1 + i;
				V(i, j) = -1.0;
				T(i, j) = TR_SS;
				W(i, j) = 0.0;
				for(int k = i; k <= j; k ++)
					W(i, j) += gamma*ssp[k];
			}
		}
		else
		{
			for(int i = 0; i <= length - l; i ++)
			{
				int j = l - 1 + i;
				if(!IsBasePair(enc[i], enc[j]))
					V(i, j) = -1.0;
				else
					V(i, j) = 2.0*(1.0 - gamma)*bppm(i, j) + W(i + 1, j - 1);
				W(i, j) = V(i, j);
				T(i, j) = TR_PAIR;
				for(int k = i; k < j; k ++)
				{
					double value = W(i, k) + W(k + 1, j);
					if(value > W(i, j))
					{
						W(i, j) = value;
						T(i, j) = k;
					}
				}
			}
		}
	}
}

void RMEPost::TraceBack(const Matrix<int>& T)
{
	stack<BasePair> S;
	string s(length, 0);
	S.push(BasePair(0, length - 1));
	while(!S.empty())
	{
		BasePair pair = S.top();
		S.pop();
		int i = pair.i;
		int j = pair.j;
		
		if(j - i > 0)
		{
			if(T(i, j) == TR_PAIR)
			{
				s[i] = '(';
				s[j] = ')';
				ct->basepr[1][i+1] = j+1;
				ct->basepr[1][j+1] = i+1;
				S.push(BasePair(i + 1, j - 1));
				
			}
			else if(T(i, j) == TR_SS)
			{
				s[i] = '.';
				s[j] = '.';
				ct->basepr[1][i+1] = 0;
				ct->basepr[1][j+1] = 0;
				S.push(BasePair(i + 1, j - 1));
			}
			else
			{
				int k = T(i, j);
				S.push(BasePair(i, k));
				S.push(BasePair(k + 1, j));
			}
		}
		//loop regions
		else if(j - i == 0)
		{
			s[i] = '.';
			ct->basepr[1][i + 1] = 0;
		}
	}
	//cout << s << endl;
}

void RMEPost::SaveAsDotBracket(const std::string& fileName)
{
	if(ct->numofstructures > 0)
		ct->writedotbracket(fileName.c_str());
	else
		Die("Cannot save as dot bracket: no structures have been predicted");
}

void RMEPost::SaveAsCt(const std::string& fileName)
{
	if(ct->numofstructures > 0)
		ctout(ct, fileName.c_str());
	else
		Die("Cannot save as dot bracket: no structures have been predicted");
}

void RMEPost::AddOptions(ParseCommandLine* parser)
{
	gamma1Options.push_back("--gamma1");
	parser->addOptionFlagsWithParameters(gamma1Options, "Parameter gamma1 for base pairing probability. Default is 0.05.");
	
	gamma2Options.push_back("--gamma2");
	parser->addOptionFlagsWithParameters(gamma2Options, "Parameter gamma2 for single-stranded probability. Default is 0.05.");
	
	gammaOptions.push_back("--gamma");
	parser->addOptionFlagsWithParameters(gammaOptions, "Parameter gamma for MaxExpect. Default is 1.0.");
	
	minHelixLenOptions.push_back("--min-helix-len");
	parser->addOptionFlagsWithParameters(minHelixLenOptions, "Minimal helix length. Default is 0.");
	
	maxHelixLenOptions.push_back("--max-helix-len");
	parser->addOptionFlagsWithParameters(maxHelixLenOptions, "Maximal helix length. Default is 4.");
}

bool RMEPost::GetOptions(ParseCommandLine* parser)
{
	if(parser->isError())
		return false;

	if(parser->contains(gamma1Options))
		parser->setOptionDouble(gamma1Options, gamma1);
	if(parser->contains(gamma2Options))
		parser->setOptionDouble(gamma2Options, gamma2);
	if(parser->contains(gammaOptions))
		parser->setOptionDouble(gammaOptions, gamma);
	if(parser->contains(minHelixLenOptions))
		parser->setOptionInteger(minHelixLenOptions, minHelixLen);
	if(parser->contains(maxHelixLenOptions))
		parser->setOptionInteger(maxHelixLenOptions, maxHelixLen);
	return true;
}