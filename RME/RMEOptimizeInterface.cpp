#include <cstdlib>
#include <sstream>
#include <iterator>
#include <map>
#include <iomanip>
using namespace std;
// multi-threading
#include <pthread.h>

// from RNAstructure package
#include <src/ErrorChecker.h>
#include <src/ParseCommandLine.h>
#include <RNA_class/RNA.h>

#include "RMEPre.h"
#include "RMEPost.h"
#include "Utils.h"


vector<double> ExpandRange(const string& s)
{
    vector<string> tokens = SplitString(s, ':');
    if(tokens.size() != 3)
        return vector<double>();

    double start = strtod(tokens[0].c_str(), NULL);
    double stop = strtod(tokens[1].c_str(), NULL);
    double step = strtod(tokens[2].c_str(), NULL);

    if((stop < start) || (step < 0))
        return vector<double>();

    vector<double> result;
    stop += step*0.5;
    for(double a = start; a <= stop; a += step)
        result.push_back(a);
    return result;
}

struct RMEOptimizeParam
{
    // files
    string inFile;
    string paramFile;
    // m
    vector<double> mList;
    double mBest;
    vector<string> mRangeOptions;
    vector<string> mListOptions;
    // gamma1
    vector<double> gamma1List;
    double gamma1Best;
    vector<string> gamma1RangeOptions;
    vector<string> gamma1ListOptions;
    // gamma2
    vector<double> gamma2List;
    double gamma2Best;
    vector<string> gamma2RangeOptions;
    vector<string> gamma2ListOptions;
    // threads
    int nThreads;
    vector<string> nThreadsOptions;
    // detail file
    bool detailFileGiven;
    string detailFile;
    vector<string> detailFileOptions;
    // verbose flag
    bool verbose;
    vector<string> verboseOptions;
    // constructor
    RMEOptimizeParam();
};

RMEOptimizeParam::RMEOptimizeParam()
{
    // default parameters
    mList = ExpandRange("0:1:0.05");
    gamma1List = ExpandRange("0:1:0.05");
    gamma2List = ExpandRange("0:1:0.05");
    nThreads = 4;

    mBest = 0.6;
    gamma1Best = 0.05;
    gamma2Best = 0.05;

    detailFileGiven = false;
    verbose = false;
}

void AddOptions(ParseCommandLine* parser, RMEOptimizeParam& params)
{
    // parameter m
    params.mRangeOptions.push_back("--m-range");
    parser->addOptionFlagsWithParameters(params.mRangeOptions,
                                         "Range of m to be scanned. Specified in start:end:step format. Not allowed with --m-list. "
                                         "Default is 0:1:0.05.");
    params.mListOptions.push_back("--m-list");
    parser->addOptionFlagsWithParameters(params.mListOptions,
                                         "List of m values delimited by \",\". Not allowed with --m-range.");

    // parameter gamma1
    params.gamma1RangeOptions.push_back("--gamma1-range");
    parser->addOptionFlagsWithParameters(params.gamma1RangeOptions,
                                         "Range of gamma1 to be scanned. Specified in start:end:step format. Not allowed with --gamma1-list. "
                                         "Default is 0:1:0.05.");
    params.gamma1ListOptions.push_back("--gamma1-list");
    parser->addOptionFlagsWithParameters(params.gamma1ListOptions,
                                         "List of gamma1 values delimited by \",\". Not allowed with --gamma1-range.");

    // parameter gamma2
    params.gamma2RangeOptions.push_back("--gamma2-range");
    parser->addOptionFlagsWithParameters(params.gamma2RangeOptions,
                                         "Range of gamma2 to be scanned. Specified in start:end:step format. Not allowed with --gamma2-list. "
                                         "Default is 0:1:0.05.");
    params.gamma2ListOptions.push_back("--gamma2-list");
    parser->addOptionFlagsWithParameters(params.gamma2ListOptions,
                                         "List of gamma2 values delimited by \",\". Not allowed with --gamma2-range.");

    // threads
    params.nThreadsOptions.push_back("-p");
    params.nThreadsOptions.push_back("--threads");
    parser->addOptionFlagsWithParameters(params.nThreadsOptions,
                                         "Number of threads to run in parallel. Default is 4.");

    // details file
    params.detailFileOptions.push_back("--details");
    parser->addOptionFlagsWithParameters(params.detailFileOptions,
                                         "Name of the file to which the details of the optimization are written. "
                                         "Defailt is not to write details.");

    params.verboseOptions.push_back("--verbose");
    parser->addOptionFlagsNoParameters(params.verboseOptions,
                                       "Verbose output. Default is off. ");
}

bool GetOptions(ParseCommandLine* parser, RMEOptimizeParam& params)
{
    // parameter m
    if(!parser->isError())
    {
        bool mRangeOptionsGiven = parser->contains(params.mRangeOptions);
        bool mListOptionsGiven = parser->contains(params.mListOptions);
        if(mRangeOptionsGiven && mListOptionsGiven)
            parser->setErrorSpecialized("--m-range is not allowed with --m-list");
        else
        {
            if(mRangeOptionsGiven)
            {
                string s = parser->getOptionString(params.mRangeOptions, false);
                params.mList = ExpandRange(s);
                if(params.mList.size() < 1)
                    parser->setError("range specification after --m-range");
            }
            if(mListOptionsGiven)
            {
                string s = parser->getOptionString(params.mListOptions, false);
                vector<string> tokens = SplitString(s, ',');
                if(tokens.size() < 1)
                    parser->setError("list specification after --m-list");
                else
                {
                    // convert string to double
                    for(size_t i = 0; i < tokens.size(); i ++)
                        params.mList.push_back(strtod(tokens[i].c_str(), NULL));
                }
            }
        }
    }

    // parameter gamma1
    if(!parser->isError())
    {
        bool gamma1RangeOptionsGiven = parser->contains(params.gamma1RangeOptions);
        bool gamma1ListOptionsGiven = parser->contains(params.gamma1ListOptions);
        if(gamma1RangeOptionsGiven && gamma1ListOptionsGiven)
            parser->setErrorSpecialized("--gamma1-range is not allowed with --gamma1-list");
        else
        {
            if(gamma1RangeOptionsGiven)
            {
                string s = parser->getOptionString(params.gamma1RangeOptions, false);
                params.gamma1List = ExpandRange(s);
                if(params.gamma1List.size() < 1)
                    parser->setError("range specification after --gamma1-range");
            }
            if(gamma1ListOptionsGiven)
            {
                string s = parser->getOptionString(params.gamma1ListOptions, false);
                vector<string> tokens = SplitString(s, ',');
                if(tokens.size() < 1)
                    parser->setError("list specification after --gamma1-list");
                else
                {
                    // convert string to double
                    for(size_t i = 0; i < tokens.size(); i ++)
                        params.gamma1List.push_back(strtod(tokens[i].c_str(), NULL));
                }
            }
        }
    }

    // parameter gamma2
    if(!parser->isError())
    {
        bool gamma2RangeOptionsGiven = parser->contains(params.gamma2RangeOptions);
        bool gamma2ListOptionsGiven = parser->contains(params.gamma2ListOptions);
        if(gamma2RangeOptionsGiven && gamma2ListOptionsGiven)
            parser->setErrorSpecialized("--gamma2-range is not allowed with --gamma2-list");
        else
        {
            if(gamma2RangeOptionsGiven)
            {
                string s = parser->getOptionString(params.gamma2RangeOptions, false);
                params.gamma2List = ExpandRange(s);
                if(params.gamma2List.size() < 1)
                    parser->setError("range specification after --gamma1-range");
            }
            if(gamma2ListOptionsGiven)
            {
                string s = parser->getOptionString(params.gamma2ListOptions, false);
                vector<string> tokens = SplitString(s, ',');
                if(tokens.size() < 1)
                    parser->setError("list specification after --gamma2-list");
                else
                {
                    // convert string to double
                    for(size_t i = 0; i < tokens.size(); i ++)
                        params.gamma2List.push_back(strtod(tokens[i].c_str(), NULL));
                }
            }
        }
    }

    // threads
    if(!parser->isError())
    {
        if(parser->contains(params.nThreadsOptions))
        {
            parser->setOptionInteger(params.nThreadsOptions, params.nThreads);
            if(params.nThreads < 1)
                Die("Invalid number of jobs: %d", params.nThreads);
        }
    }

    if(!parser->isError())
    {
        if(parser->contains(params.detailFileOptions))
        {
            params.detailFile = parser->getOptionString(params.detailFileOptions, false);
            params.detailFileGiven = true;
        }
    }

    if(!parser->isError())
    {
        if(parser->contains(params.verboseOptions))
            params.verbose = true;
    }

    params.inFile = parser->getParameter(1);
    params.paramFile = parser->getParameter(2);

    return !parser->isError();
}

std::vector<structure*> CreateRefCtList(std::vector<RMEInputDataRecord>& inputData)
{
    vector<structure*> refCtList(inputData.size());
    for(size_t i = 0; i < inputData.size(); i ++)
    {
        // construct a ct object
        structure* ct = new structure(2);
        ct->numofbases = inputData[i].length;
        ct->allocate(inputData[i].length);
        ct->numofstructures = 1;
        ct->checknumberofstructures();
        ct->energy[1] = 0;
        for(int j = 1; j <= inputData[i].length; j ++)
        {
            ct->nucs[j] = inputData[i].seq[j - 1];
            ct->basepr[1][j] = inputData[i].refbp[j - 1];
        }
        strcpy(ct->ctlabel[1], inputData[i].name.c_str());
        refCtList[i] = ct;
    }
    return refCtList;
}

void FreeRefCtList(std::vector<structure*>& refCtList)
{
    for(size_t i = 0; i < refCtList.size(); i ++)
        delete refCtList[i];
}

struct RMEOptimizeItem
{
    // input
    RNA* strand;
    structure* refct;
    TProgressBar<int>* progress;

    string name;
    string seq;
    list<RMEProbRecord> probRecords;

    double m;
    double gamma1;
    double gamma2;


    // output
    double sensitivity;
    double ppv;
    double arithmean;

    bool available;
    RMEOptimizeItem(): strand(NULL), refct(NULL), progress(NULL),
        available(false) {}
};


// input: seq, probRecords, m, gamma1, gamma2
// output: sensitivity, ppv
void RunRME(void* arg)
{
    RMEOptimizeItem* item = (RMEOptimizeItem*)arg;
    // run RME-Pre
    RMEPre pre(item->seq, item->probRecords);
    pre.m = item->m;
    if(!pre.CalcPartitionFunction())
        Die("Failed to calculate partition function");
    RNA* strand = pre.GetStrand();
    // run RME-Post
    RMEPost post(strand, item->probRecords);
    post.gamma1 = item->gamma1;
    post.gamma2 = item->gamma2;
    post.Fold();
    structure* predct = post.GetStructure();
    // get score
    GetScore(predct, item->refct, item->sensitivity, item->ppv);
    item->arithmean = (item->sensitivity + item->ppv)/2.0;
    item->available = true;

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

// input: strand, probRecords, m, gamma1, gamma2
// output: strand
void RunRMEPreWithStrand(void* arg)
{

    RMEOptimizeItem* item = (RMEOptimizeItem*)arg;
    // run RME-Pre
    RMEPre pre(item->strand, item->probRecords);
    pre.m = item->m;
    if(!pre.CalcPartitionFunction())
        Die("Failed to calculate partition function");
    item->strand = pre.GetStrand();

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

// input: strand, probRecords, m, gamma1, gamma2
// output: ppv, sensitivity
void RunRMEPostWithStrand(void* arg)
{
    RMEOptimizeItem* item = (RMEOptimizeItem*)arg;
    // run RME-Post
    RMEPost post(item->strand, item->probRecords);
    post.gamma1 = item->gamma1;
    post.gamma2 = item->gamma2;
    post.Fold();
    structure* predct = post.GetStructure();
    // get score
    GetScore(predct, item->refct, item->sensitivity, item->ppv);
    item->arithmean = (item->sensitivity + item->ppv)/2.0;

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

class DetailWriter
{
public:
    DetailWriter(const string& fileName): _out(fileName.c_str())
    {
        if(!_out)
            Die("Cannot write to detail file %s", fileName.c_str());
    }
    void Write(const std::vector<RMEOptimizeItem>& items, int nSeqs, int nParams);
    void WriteTitle(const std::string& title);
    void Close() {
        _out.close();
    }
private:
    std::ofstream _out;
};

void DetailWriter::Write(const std::vector<RMEOptimizeItem>& items, int nSeqs, int nParams)
{
    _out << "m\tgamma1\tgamma2\tSensitivity\tPPV\tArithMean" << endl;
    for(int iParam = 0; iParam < nParams; iParam ++)
    {
        double avgscore = 0.0;
        double avgsensitivity = 0.0;
        double avgppv = 0.0;
        for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
        {
            int index = iParam*nSeqs + iSeq;
            avgscore += items[index].arithmean;
            avgsensitivity += items[index].sensitivity;
            avgppv += items[index].ppv;
            /*
            _out << items[index].name << "\t"
            	<< items[index].m << "\t"
            	<< items[index].gamma1 << "\t"
            	<< items[index].gamma2 << "\t"
            	<< items[index].sensitivity << "\t"
            	<< items[index].ppv << "\t"
            	<< items[index].arithmean << endl;*/
        }
        avgscore /= nSeqs;
        avgsensitivity /= nSeqs;
        avgppv /= nSeqs;

        _out << items[iParam*nSeqs].m << "\t"
             << items[iParam*nSeqs].gamma1 << "\t"
             << items[iParam*nSeqs].gamma2 << "\t"
             << avgsensitivity << "\t"
             << avgppv << "\t"
             << avgscore << endl;
    }
}

void DetailWriter::WriteTitle(const std::string& title)
{
    _out << "=== " << title << " ===" << endl;
}

void OptimizeRMEPreMT(std::vector<RMEInputDataRecord>& inputData,
                      std::vector<structure*>& refCtList,
                      RMEOptimizeParam& params,
                      DetailWriter* detailWriter)
{
    cout << "Optimizing parameters for RME-Pre" << endl;

    int nParams = params.mList.size();
    int nSeqs = inputData.size();

    // progress bar
    TProgressBar<int> progress(nSeqs*nParams, 1);
    // data for threads
    vector<RMEOptimizeItem> items(nParams*nSeqs);
    for(int iParam = 0; iParam < nParams; iParam ++)
    {
        for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
        {
            int index = iParam*nSeqs + iSeq;
            items[index].name = inputData[iSeq].name;
            items[index].m = params.mList[iParam];
            items[index].gamma1 = 0;
            items[index].gamma2 = 0;
            items[index].seq = inputData[iSeq].seq;
            items[index].probRecords = inputData[iSeq].probRecords;
            items[index].refct = refCtList[iSeq];
            items[index].progress = &progress;
        }
    }
    // run RME
    progress.SetCaption("Optimizing RME-Pre");
    progress.Show();
    ApplyFunction<RMEOptimizeItem>(RunRME, items, params.nThreads);
    cout << endl;

    // find best parameters
    double maxscore = 0.0;
    for(int iParam = 0; iParam < nParams; iParam ++)
    {
        double avgscore = 0.0;
        for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
        {
            int index = iParam*nSeqs + iSeq;
            avgscore += items[index].arithmean;
        }
        if(avgscore > maxscore)
        {
            maxscore = avgscore;
            params.mBest = items[iParam*nSeqs].m;
        }
    }

    cout << "Optimized parameters: m = " << params.mBest << endl;

    // write the optimization details to file
    if(detailWriter)
    {
        detailWriter->WriteTitle("RME-Pre optimization results");
        detailWriter->Write(items, nSeqs, nParams);
    }
}

void OptimizeRMEPostMT(std::vector<RMEInputDataRecord>& inputData,
                       std::vector<structure*>& refCtList,
                       RMEOptimizeParam& params,
                       DetailWriter* detailWriter)
{
    cout << "Optimizing parameters for RME-Post" << endl;

    int nParams = params.gamma1List.size()*params.gamma2List.size();
    int nSeqs = inputData.size();

    // progress bar
    TProgressBar<int> progress(nSeqs, 1);
    // allocate strands
    vector<RMEOptimizeItem> preItems(nSeqs);
    for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
    {
        preItems[iSeq].m = params.mBest;
        preItems[iSeq].strand = new RNA(inputData[iSeq].seq.c_str(), true);
        preItems[iSeq].probRecords = inputData[iSeq].probRecords;
        // set progress
        preItems[iSeq].progress = &progress;
    }
    // Run RME-Pre with best m
    cout << "Run RME-Pre with best parameters: m = " << params.mBest << endl;
    progress.SetCaption("Running RME-Pre");
    progress.Show();
    ApplyFunction<RMEOptimizeItem>(RunRMEPreWithStrand, preItems, params.nThreads);
    cout << endl;

    progress.Reset(nSeqs*nParams, 1);
    // data for RME-Post
    vector<RMEOptimizeItem> items(nParams*nSeqs);
    for(int iParam = 0; iParam < nParams; iParam ++)
    {
        for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
        {
            int index = iParam*nSeqs + iSeq;
            items[index].name = inputData[iSeq].name;
            items[index].m = params.mBest;
            items[index].gamma1 = params.gamma1List[iParam/params.gamma2List.size()];
            items[index].gamma2 = params.gamma2List[iParam%params.gamma2List.size()];
            items[index].seq = inputData[iSeq].seq;
            items[index].probRecords = inputData[iSeq].probRecords;
            items[index].refct = refCtList[iSeq];
            items[index].strand = preItems[iSeq].strand;
            // set progress
            items[index].progress = &progress;
        }
    }
    // Run RME-Post
    progress.SetCaption("Optimizing RME-Post");
    progress.Show();
    ApplyFunction<RMEOptimizeItem>(RunRMEPostWithStrand, items, params.nThreads);
    cout << endl;
    // free strands
    for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
        delete preItems[iSeq].strand;
    // find best parameters
    double maxscore = 0.0;
    for(int iParam = 0; iParam < nParams; iParam ++)
    {
        double avgscore = 0.0;
        for(int iSeq = 0; iSeq < nSeqs; iSeq ++)
        {
            int index = iParam*nSeqs + iSeq;
            avgscore += items[index].arithmean;
        }
        if(avgscore > maxscore)
        {
            maxscore = avgscore;
            params.gamma1Best = items[iParam*nSeqs].gamma1;
            params.gamma2Best = items[iParam*nSeqs].gamma2;
        }
    }

    cout << "Optimized parameters: gamma1 = " << params.gamma1Best
         << ", gamma2 = " << params.gamma2Best << endl;

    // write the optimization details to file
    if(detailWriter)
    {
        detailWriter->WriteTitle("RME-Post optimization results");
        detailWriter->Write(items, nSeqs, nParams);
    }
}

int main(int argc, char** argv)
{
    // Create the command line parser and build in its required parameters.
    ParseCommandLine* parser = new ParseCommandLine("RME-Optimize");

    // required parameters
    parser->addParameterDescription("input file",
                                    "The name of the file containing RNA sequence. "
                                    "constraint base pairing probability and reference structure. "
                                    "This file should include 5 columns and each line specifies a base. "
                                    "Column 1: RNA sequence name. "
                                    "Column 2: Position of the base in the RNA sequence (1-based). "
                                    "Column 3: Base-pairing probability transformed from structural probing data. "
                                    "Column 4: Name of the base (A, C, G, U). "
                                    "Column 5: Reference structure. 0 for single-stranded. Non-zero is the position of the base to which it pairs." );
    parser->addParameterDescription("parameter file",
                                    "The name of the file to which optimized parameters will be written.");

    RMEOptimizeParam params;
    AddOptions(parser, params);
    // parse command line
    parser->parseLine(argc, argv);

    if(GetOptions(parser, params))
    {
        if(!getenv("DATAPATH"))
            Die("environment variable $DATAPATH should be set to {RNAstructure root}/data_tables/");

        if(params.verbose)
        {
            cout << "m list: ";
            for(size_t i = 0; i < params.mList.size(); i ++)
                cout << params.mList[i] << " ";
            cout << endl;
            cout << "gamma1 list: ";
            for(size_t i = 0; i < params.gamma1List.size(); i ++)
                cout << params.gamma1List[i] << " ";
            cout << endl;
            cout << "gamma2 list: ";
            for(size_t i = 0; i < params.gamma2List.size(); i ++)
                cout << params.gamma2List[i] << " ";
            cout << endl;
        }
        // read input data
        vector<RMEInputDataRecord> inputData = ReadRMEInputData(params.inFile, true);

        if(params.verbose)
        {
            for(size_t i = 0; i < inputData.size(); i ++)
            {
                cout << inputData[i].name << ", length: " << inputData[i].length;
                cout << ", number of positions that have probability data: " << inputData[i].probRecords.size() << endl;
                //cout << inputData[i].seq << endl;
            }
        }

        // open parameter file
        ofstream foutParam(params.paramFile.c_str());
        if(!foutParam)
            Die("Cannot write to parameter file %s", params.paramFile.c_str());

        // open detail file
        DetailWriter* detailWriter = NULL;
        if(params.detailFileGiven)
            detailWriter = new DetailWriter(params.detailFile);

        vector<structure*> refCtList = CreateRefCtList(inputData);
        cout << "Multi-threading enabled, using " << params.nThreads << " threads." << endl;
        OptimizeRMEPreMT(inputData, refCtList, params, detailWriter);
        OptimizeRMEPostMT(inputData, refCtList, params, detailWriter);

        // write parameters to file
        cout << "Writing optimized parameters to file: " << params.paramFile << endl;
        foutParam << "m\t" << params.mBest << endl;
        foutParam << "gamma1\t" << params.gamma1Best << endl;
        foutParam << "gamma2\t" << params.gamma2Best << endl;
        foutParam.close();
        // detail file
        if(detailWriter)
            detailWriter->Close();

        FreeRefCtList(refCtList);
    }

    return 0;
}