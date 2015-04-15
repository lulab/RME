// AveragePairProb.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cstdlib>
#include <list>
#include <string>
using namespace std;

#include <RNA_class/RNA.h>

#include "Utils.h"

struct AveragePairProbItem
{
    string seqfile;
    double result;
};

double AveragePairProb(const char* seqfile)
{
    //Correct # of parameters, try the calculation
    //Allocate RNA, 2 = sequence file as input
    RNA* rna = new RNA(seqfile, 2);
    double result = -1;
    //Check the error codes to make sure sequence was read correctly
    if (rna->GetErrorCode()==0)
    {
        //Now call the partition function
        int error = rna->PartitionFunction();
        if (error==0) {
        //No problem with partition function calculation
        vector<double> pairing(rna->GetSequenceLength()+1,0.0);
        double mean=0;
        //Note that sequences are 1 indexed in the RNA class
        for (int i=1; i < rna->GetSequenceLength(); ++i)
        {
            for (int j=i+1; j <= rna->GetSequenceLength(); ++j)
            {
                //Add up the pair probabilities for each nucleotide, i to the 3' direction
                pairing[i]+=rna->GetPairProbability(i,j);
            }//end for int j...
            //NOTE: I calculated the individual pair probabilities, although that was not really necessary.
            //I did the long calculation to make clear how the per nucleotide pair probabilities could be calculated.
            mean+=pairing[i];
        }//end for (int i...
            mean = 2 * mean / ((double) (rna->GetSequenceLength()));
            result = mean;
        }
        else
        {
            //There was an error in the parition function calculation, report it
            cerr << rna->GetErrorMessage(error);
        }
    }//end rna->GetErrorCode()==0
    else
    {
        //There was an error when reading the sequence, report the error
        cerr << rna->GetErrorMessage(rna->GetErrorCode());
    }
    //Clean up memory
    delete rna;
    
    return result;
}

void RunAveragePairProb(void* args)
{
    AveragePairProbItem* item = (AveragePairProbItem*)args;
    item->result = AveragePairProb(item->seqfile.c_str());
}

void RunAveragePairProbMT(const vector<string>& seqfiles, int nThreads)
{
    vector<AveragePairProbItem> items(seqfiles.size());
    for(int i = 0; i < items.size(); i ++)
        items[i].seqfile = seqfiles[i];
    
    ApplyFunction<AveragePairProbItem>(RunAveragePairProb, items, nThreads);
    
    // output results
    for(int i = 0; i < items.size(); i ++)
        cout << items[i].seqfile << "\t" << items[i].result << endl;
}

void mainMT(int argc, char** argv)
{
    if (argc < 2) 
    {
        //Provide usage information
        cerr << "USAGE: AveragePairProb [-p threads] seqfile1 [seqfile2 ...]" << endl;
        exit(-1);
    }
    else
    {
        int argi = 1;
        int nThreads = 1;
        if(!strcmp(argv[argi], "-p"))
        {
            if(argc < 3)
            {
                cerr << "Error: -p: missing an argument" << endl;
                exit(-1);
            }
            argi ++;
            nThreads = atoi(argv[argi]);
            if(nThreads <= 0)
            {
                cerr << "Error: -p: invalid argument value " << argv[argi] << endl;
                exit(-1);
            }
            argi ++;
        }
        if(argi >= argc)
        {
            cerr << "Error: missing sequence file names" << endl;
            exit(-1);
        }
        
        vector<string> seqfiles;
        while(argi < argc)
        {
            if(argv[argi][0] == '-')
            {
                cerr << "Error: " << argv[argi] << ": unrecognized option" << endl;
                exit(-1);
            }
            seqfiles.push_back(string(argv[argi]));
            argi ++;
        }
        for(int i = 0; i < seqfiles.size(); i ++)
        {
            if(!FileExists(seqfiles[i]))
            {
                cerr << "Error: file does not exist: " << seqfiles[i] << endl;
                exit(-1);
            }
        }
        RunAveragePairProbMT(seqfiles, nThreads);
    }
}

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        cerr << "Usage: " << argv[0] << " seqfile" << endl;
        return -1;
    }
    double result = AveragePairProb(argv[1]);
    cout << result << endl;
    
    return 0;
}

