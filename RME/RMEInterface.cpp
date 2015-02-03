#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
using namespace std;
#include <string.h>
#include <errno.h>
//from RNAstructure
#include "../RNA_class/RNA.h"
#include "../src/ParseCommandLine.h"

#include "Utils.h"
#include "RMEPre.h"
#include "RMEPost.h"

struct RMEParam
{
    string inFile;
    string outDir;

    // parameter file
    std::string paramFile;
    std::vector<std::string> paramFileOptions;
    // threads
    int nThreads;
    std::vector<std::string> nThreadsOptions;
    // epsilon
    double epsilon;
    std::vector<std::string> epsilonOptions;
    // m
    double m;
    std::vector<std::string> mOptions;
    // temperature
    double temperature;
    std::vector<std::string> tempOptions;
    // gamma1
    double gamma1;
    std::vector<std::string> gamma1Options;
    // gamma2
    double gamma2;
    std::vector<std::string> gamma2Options;
    // gamma
    double gamma;
    std::vector<std::string> gammaOptions;
    // min helix length
    int minHelixLen;
    std::vector<std::string> minHelixLenOptions;
    // max helix length
    int maxHelixLen;
    std::vector<std::string> maxHelixLenOptions;
    // initialization
    RMEParam():
        nThreads(4),
        epsilon(0.01),
        m(0.5),
        temperature(310.25),
        gamma1(0.1),
        gamma2(0.1),
        gamma(0.5),
        minHelixLen(0),
        maxHelixLen(4) {}
};

void AddOptions(ParseCommandLine* parser, RMEParam& params)
{
    params.epsilonOptions.push_back("--epsilon");
    parser->addOptionFlagsWithParameters(params.epsilonOptions, "A small value used when the constraint probability falls below 0. "
                                         "Default is 0.01.");

    params.mOptions.push_back("-m");
    params.mOptions.push_back("--m");
    parser->addOptionFlagsWithParameters(params.mOptions, "Constraint weight. Default is 0.5.");

    params.tempOptions.push_back("-t");
    params.tempOptions.push_back("-T");
    params.tempOptions.push_back("--temperature" );
    parser->addOptionFlagsWithParameters(params.tempOptions, "Specify the temperature at which calculation takes place in Kelvin. "
                                         "Default is 310.15 K, which is 37 degrees C." );

    params.gamma1Options.push_back("--gamma1");
    parser->addOptionFlagsWithParameters(params.gamma1Options, "Parameter gamma1 for base pairing probability. Default is 0.1.");

    params.gamma2Options.push_back("--gamma2");
    parser->addOptionFlagsWithParameters(params.gamma2Options, "Parameter gamma2 for single-stranded probability. Default is 0.1.");

    params.gammaOptions.push_back("--gamma");
    parser->addOptionFlagsWithParameters(params.gammaOptions, "Parameter gamma for MaxExpect. Default is 0.5.");

    params.minHelixLenOptions.push_back("--min-helix-len");
    parser->addOptionFlagsWithParameters(params.minHelixLenOptions, "Minimal helix length. Default is 0.");

    params.maxHelixLenOptions.push_back("--max-helix-len");
    parser->addOptionFlagsWithParameters(params.maxHelixLenOptions, "Maximal helix length. Default is 4.");

    params.nThreadsOptions.push_back("-p");
    params.nThreadsOptions.push_back("--threads");
    parser->addOptionFlagsWithParameters(params.nThreadsOptions,
                                         "Number of threads to run in parallel. Default is 4.");

    params.paramFileOptions.push_back("-i");
    params.paramFileOptions.push_back("--params");
    parser->addOptionFlagsWithParameters(params.paramFileOptions,
                                         "Name of the file containing parameters. Each line contains a (name, value) pair. ");
}

void ReadParamFile(std::string& fileName, RMEParam& params)
{
    ifstream fin(fileName.c_str());
    while(!fin.eof())
    {
        string line;
        std::getline(fin, line);

        string name;
        double value;
        istringstream is(line);
        is >> name >> value;

        if(name == "m")
            params.m = value;
        else if(name == "gamma1")
            params.gamma1 = value;
        else if(name == "gamma2")
            params.gamma2 = value;
        else if(name == "gamma")
            params.gamma = value;
    }
    fin.close();
}

bool GetOptions(ParseCommandLine* parser, RMEParam& params)
{
    if(!parser->isError() && parser->contains(params.paramFileOptions))
    {
        params.paramFile = parser->getOptionString(params.paramFileOptions, true);
        ReadParamFile(params.paramFile, params);
    }

    if(!parser->isError() && parser->contains(params.tempOptions))
    {
        parser->setOptionDouble(params.tempOptions, params.temperature);
        if(params.temperature < 0)
            parser->setError("temperature");
    }

    if(!parser->isError() && parser->contains(params.epsilonOptions))
    {
        parser->setOptionDouble(params.epsilonOptions, params.epsilon);
        if(params.epsilon <= 0)
            parser->setError("epsilon");
    }

    if(!parser->isError() && parser->contains(params.mOptions))
    {
        parser->setOptionDouble(params.mOptions, params.m);
        if(params.m <= 0)
            parser->setError("m");
    }

    if(!parser->isError() && parser->contains(params.gamma1Options))
    {
        parser->setOptionDouble(params.gamma1Options, params.gamma1);
        if(params.gamma1 < 0 || params.gamma1 > 1)
            parser->setError("gamma1");
    }

    if(!parser->isError() && parser->contains(params.gamma2Options))
    {
        parser->setOptionDouble(params.gamma2Options, params.gamma2);
        if(params.gamma2 < 0 || params.gamma2 > 1)
            parser->setError("gamma2");
    }

    if(!parser->isError() && parser->contains(params.gammaOptions))
    {
        parser->setOptionDouble(params.gammaOptions, params.gamma);
        if(params.gamma < 0 || params.gamma > 1)
            parser->setError("gamma");
    }

    if(!parser->isError() && parser->contains(params.minHelixLenOptions))
    {
        parser->setOptionInteger(params.minHelixLenOptions, params.minHelixLen);
        if(params.minHelixLen < 0)
            parser->setError("minimum helix length");
    }

    if(!parser->isError() && parser->contains(params.maxHelixLenOptions))
    {
        parser->setOptionInteger(params.maxHelixLenOptions, params.maxHelixLen);
        if(params.maxHelixLen < 0)
            parser->setError("maximum helix length");
    }

    if(!parser->isError() && parser->contains(params.nThreadsOptions))
    {
        parser->setOptionInteger(params.nThreadsOptions, params.nThreads);
        if(params.nThreads < 1)
            parser->setError("number of threads");
    }

    return !parser->isError();
}

struct RMEItem
{
    // input
    RNA* strand;
    TProgressBar<int>* progress;

    string name;
    string seq;
    string ctFile;
    list<RMEProbRecord> probRecords;
    RMEParam params;

    RMEItem(): strand(NULL), progress(NULL)
    {}
};

void RunRME(void* arg)
{
    RMEItem* item = (RMEItem*)arg;
    // Run RME-Pre
    RMEPre pre(item->seq, item->probRecords);
    pre.m = item->params.m;
    pre.temperature = item->params.temperature;
    pre.epsilon = item->params.epsilon;
    if(!pre.CalcPartitionFunction())
    {
        cerr << "Failed to calculate partition function for sequence: " << item->name << endl;
        return;
    }
    item->strand = pre.GetStrand();
    // Run RME-Post
    RMEPost post(item->strand, item->probRecords);
    post.SetName(item->name);
    post.gamma1 = item->params.gamma1;
    post.gamma2 = item->params.gamma2;
    post.gamma = item->params.gamma;
    post.minHelixLen = item->params.minHelixLen;
    post.maxHelixLen = item->params.maxHelixLen;
    post.Fold();
    post.SaveAsCt(item->ctFile);
    // show progress
    if(item->progress)
    {
        static pthread_mutex_t mutex_prog = PTHREAD_MUTEX_INITIALIZER;
        pthread_mutex_lock(&mutex_prog);
        item->progress->Increment();
        item->progress->Show();
        pthread_mutex_unlock(&mutex_prog);
    }
}

void RunRMEMT(std::vector<RMEInputDataRecord>& inputData,
              RMEParam& params)
{
    int nSeqs = inputData.size();
    vector<RMEItem> items(nSeqs);
    // progress bar
    TProgressBar<int> progress(nSeqs, 1);

    for(int i = 0; i < nSeqs; i ++)
    {
        items[i].name = inputData[i].name;
        items[i].ctFile = params.outDir + "/" + items[i].name + ".ct";
        items[i].seq = inputData[i].seq;
        items[i].probRecords = inputData[i].probRecords;
        items[i].params = params;
        items[i].progress = &progress;
    }
    // run RME
    progress.SetCaption("Running RME");
    progress.Show();
    ApplyFunction<RMEItem>(RunRME, items, params.nThreads);
    cout << endl;
}

int main(int argc, char** argv)
{
    RMEParam params;
    // Create the command line parser and build in its required parameters.
    ParseCommandLine* parser = new ParseCommandLine("RME");
    // required arguments
    parser->addParameterDescription("input file", "The name of a file containing input data."
                                    "This file should include 4 columns and each line specifies a base. "
                                    "Column 1: RNA sequence name. "
                                    "Column 2: Position of the base in the RNA sequence (1-based). "
                                    "Column 3: Base-pairing probability transformed from structural probing data. "
                                    "Column 4: Name of the base (A, C, G, U). ");
    parser->addParameterDescription("output directory", "The name of the directory where .ct files will be written. ");

    AddOptions(parser, params);
    // parser command line
    parser->parseLine(argc, argv);
    // get required parameters from command line arguments
    //string seqFile = parser->getParameter(1);
    //string probFile = parser->getParameter(2);
    //string ctFile = parser->getParameter(3);
    params.inFile = parser->getParameter(1);
    params.outDir = parser->getParameter(2);
    // get options
    if(GetOptions(parser, params))
    {
        // read input data without reference structures
        vector<RMEInputDataRecord> inputData = ReadRMEInputData(params.inFile, false);
        // create output directory
        if(!DirExists(params.outDir))
        {
            if(!MakeDir(params.outDir))
                Die("Cannot create output directory: %s", strerror(errno));
        }
        RunRMEMT(inputData, params);
    }
    return 0;
}