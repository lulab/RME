# RME

A tool for RNA secondary structure prediction with various types of experimental probing data. It can use experimental pairing probabilities to restrain the Partition Function, and predict the structure with maximum restrained expected accuray based on a MEA algorithm, MaxExpect (Lu *et al*., 2009 *RNA* ). It is based on the *RNAstructure* package. It also provides example scripts for transforming the SHAPE, PARS and DMS-seq data into pairing probability according a Posterior Probabilistic Model. Moreover, it also contains a utility for optimizing the parameters of RME by RME-Optimize.

For updates, please refer to https://github.com/lulab/RME


## Citation

Yang Wu, Binbin Shi *et al*., (2015) Improved prediction of RNA secondary structure by integrating the free energy model with restraints derived from experimental probing data.


## License

RME is free for non-commercial research. For commercial use,please contact the authors.


## Prerequisites

1. Linux
2. Perl (version >= 5.10.1)
3. R (version >= 3.1.0)
4. Bioconductor 
5. R packages: preprocessCore argparse mixtools reshape MASS evir


## Installation

First enter the directory of RME and type:
```
make
```
Then the executables will be placed in the `bin` directory.
You can add the `bin` directory to your `PATH` variable for convinient use:
```bash
export PATH=$PATH:/path/to/RME/bin
```

If you have moved the already built RME to another directory, make sure to rebuild RME by typing:
```
make clean
make
```
Optionally, you can install RME to another location by typing:
```
make install PREFIX=/path/to/install
```

## Usage of RME


### Step 1: prepare data files

#### a. Prepare reference structures for training RNAs

Example files for reference structures:
```
example/dat/structure/*.ct
```

#### b. Prepare data files for training and testing RNAs

Example files for training sets:
```
example/dat/data/SHAPE.train.23SrRNA.data
example/dat/data/PARS.train.25SrRNA.data
example/dat/data/DMSseq.train.25SrRNA.data
```
Example files for test sets:
```
example/dat/data/SHAPE.test.data
example/dat/data/PARS.test.data
example/dat/data/DMSseq.test.data
```
#### c. Format of data files

For 1-dimensional probing data (i.e. SHAPE data), training files contains 5 columns separated by tabs or spaces.
- column 1    name of RNA
- column 2    index on the RNA, 1-based
- column 3    reactivity on this base
- column 4    base information in capitals, ACGU
- column 5    structure information for this base, 0 stands for single-standed bases, 1 strands for paired bases

For 2-dimensional probing data (i.e. PARS and DMS-seq data), training files contains 6 columns.
- column 1    name of RNA
- column 2    index on the RNA, 1-based
- column 3    read count 1 for this base, V1 for PARS data and control for DMS-seq data
- column 4    read count 2 for this base, S1 for PARS data and vivo for DMS-seq data
- column 5    base information in capitals, ACGU
- column 6    structure information for this base, 0 stands for single-standed bases, 1 strands for paired bases

Files for test RNAs are similar with training files except that the last column (structure information) is missing.
Note that bases with probing signals are not included in these files. 



### Step 2: transform structure probing data into pairing probability  

#### a. Run preprocessing scripts

The transforming scripts accepts 4 parametesr: a training file and a test file prepared in step 1, directory for the structure files, as well as the output directory.
```bash
processSHAPE example/dat/data/SHAPE.train.23SrRNA.data example/dat/data/SHAPE.test.data example/dat/structure example/1.processing-data
processPARS example/dat/data/PARS.train.25SrRNA.data example/dat/data/PARS.test.data example/dat/structure example/1.processing-data
processDMSseq example/dat/data/DMSseq.train.25SrRNA.data example/dat/data/DMSseq.test.data example/dat/structure example/1.processing-data
```
#### b. Output files of step 2

Outfiles that will be used in step 3, which contains pairing probabilities and structure information of the training RNAs
```
example/1.processing-data/SHAPE.for-optimize-parameter.txt
example/1.processing-data/PARS.for-optimize-parameter.txt
example/1.processing-data/DMSseq.for-optimize-parameter.txt
```
Outfiles that will be used in step 4, which contains pairing probabilities of the test RNAs
```
example/1.processing-data/SHAPE.for-test.txt
example/1.processing-data/PARS.for-test.txt
example/1.processing-data/DMSseq.for-test.txt
```
### Step 3: optimize the paramter for RME

#### a. Parameter optimization by RME-Optimize
```bash
mkdir example/2.optimize-param/
cd example/2.optimize-param/
RME-Optimize ../1.processing-data/SHAPE.for-optimize-parameter.txt SHAPE.params.txt --m-range 0:1:0.05 --details SHAPE.details.txt -p 4
RME-Optimize ../1.processing-data/PARS.for-optimize-parameter.txt PARS.params.txt --m-range 0:1:0.02 --details PARS.details.txt -p 4
RME-Optimize ../1.processing-data/DMSseq.for-optimize-parameter.txt DMSseq.params.txt --m-range 0:1:0.02 --details DMSseq.details.txt -p 4
cd ../..
```
You can get more information about usage of `RME-Optimize` by typing `RME-Optimize -h`.

#### b. Output files of step 3

Output files containing the optimal parameters, which will be used in step 4
```
example/2.optimize-param/SHAPE.params.txt
example/2.optimize-param/PARS.params.txt
example/2.optimize-param/DMSseq.params.txt
```
Output files recording the detail of the optimizing process
```
example/2.optimize-param/SHAPE.details.txt
example/2.optimize-param/PARS.details.txt
example/2.optimize-param/DMSseq.details.txt
```
### Step 4: Structure prediction by RME

#### a. Structure prediction by RME
```bash
mkdir example/3.prediction/
cd example/3.prediction/
RME ../1.processing-data/SHAPE.for-test.txt SHAPE -i ../2.optimize-param/SHAPE.params.txt -p 4
RME ../1.processing-data/PARS.for-test.txt PARS -i ../2.optimize-param/PARS.params.txt -p 4
RME ../1.processing-data/DMSseq.for-test.txt DMSseq -i ../2.optimize-param/DMSseq.params.txt -p 4
cd ../..
```
You can get more information about usage of `RME` by typing `RME -h`.

#### b. Output files of the step 4

Outfiles of the predicted structures in CT format
```
example/3.prediction/SHAPE/*.ct
example/3.prediction/PARS/*.ct
example/3.prediction/DMSseq/*.ct
```

## Output files

structure prediction files

```
example/3.prediction/*/*.ct
```

pairing probability files

```
example/1.processing-data/*/*.prob
```

