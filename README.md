# RME

A tool for RNA secondary structure prediction with multiple types of experimental probing data. It can use experimental pairing probabilities to restrain the partition function, and predict the structure with maximum restrained expected accuray based on a MEA algorithm, MaxExpect (Lu *et al*., 2009 *RNA* ). It is based on the *RNAstructure* (http://rna.urmc.rochester.edu/RNAstructure.html) package. It also provides preprocessing scripts for transforming the SHAPE, PARS and DMS-seq data into pairing probability according a posterior probabilistic model. Moreover, it also contains a utility for optimizing the parameters of RME by RME-Optimize.

For updates, please refer to https://github.com/lulab/RME


## Citation

Yang Wu, Binbin Shi *et al*., (2015) Improved prediction of RNA secondary structure by integrating the free energy model with restraints derived from experimental probing data.


## License

RME is free for non-commercial research. For commercial use, please contact the authors.


## Download

Please download the software from https://github.com/lulab/RME/releases


## Prerequisites

1. Linux
2. Perl (version >= 5.10.1)
3. R (version >= 3.1.0)
4. R packages: [argparse](http://cran.r-project.org/web/packages/argparse/index.html) [mixtools](http://cran.r-project.org/web/packages/mixtools/index.html) [reshape](http://cran.r-project.org/web/packages/reshape/index.html) [MASS](http://cran.r-project.org/web/packages/MASS/index.html) [evir](http://cran.r-project.org/web/packages/evir/index.html)
5. Bioconductor package: [preprocessCore](http://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html)


## Installation

First of all, enter the directory of RME and type:
```bash
make
```
Then the executables will be placed in the `bin` directory.
You can add the `bin` directory to your `PATH` variable for convinient use:
```bash
export PATH=$PATH:/path/to/RME/bin
```

Optionally, you can install RME to another location by typing:
```bash
make install PREFIX=/path/to/install
```


## Usage of RME


### Step 1: prepare data files

#### a. Prepare reference structures for each RNA

Example files for reference structures:
```
example/dat/structure/*.ct
```

All structure files are in CT format (http://rna.urmc.rochester.edu/Text/File_Formats.html). And the files are nominated by the RNA name (also appeared in the data files) and a .ct suffix. For test RNAs without reference structures, the 5th column in the CT file can be meaningless (e.g. set to be a consistent number). See `example/dat/structure/README` for details.


#### b. Prepare sequence files for each RNA

Example files for RNA sequences
```
example/dat/sequence/*.fa
```
All data files are in FASTA format (http://rna.urmc.rochester.edu/Text/File_Formats.html). And the files are nominated by the RNA name (also appeared in the data files) and a .fa suffix. See `example/dat/sequence/README` for details.


#### c. Prepare data files for each RNA

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

#### d. Format of data files

For SHAPE data with one reactivity value for each nucleotide, training files contain 5 columns separated by a TAB. The first lines of the files specify the column name (given in parentheses).

| Column | Description |
|:----:| ------------- |
| 1 |  name of RNA (RNA)|
| 2 |  index on the RNA, 1-based (Index)|
| 3 |  reactivity on this base (Reactivity)|
| 4 |  base information in capitals, ACGU (Base)|
| 5 |  structure information for this base, 0 stands for single-standed bases, 1 strands for paired bases (Structure)|

For PARS and DMS-seq data with two read counts for each nucleotdie, training files contain 6 columns by a TAB. The first lines of the files specify the column name (given in parentheses).

| Column | Description |
|:----:| ------------- |
| 1 |  name of RNA (RNA)|
| 2 |  index on the RNA, 1-based (Index)|
| 3 |  read count 1 for this base, V1 for PARS data and control for DMS-seq data (V1 or Control)|
| 4 |  read count 2 for this base, S1 for PARS data and vivo for DMS-seq data (S1 or Vivo)|
| 5 |  base information in capitals, ACGU (Base)|
| 6 |  structure information for this base, 0 stands for single-standed bases, 1 strands for paired bases (Structure)|

Files for test RNAs are similar with the training files except that the last column (structure information) is missing.
Note that bases without probing signals are not included in these files. See `example/dat/data/README` for details.


### Step 2: transform structure probing data into pairing probability  

#### a. Run preprocessing scripts

The preprocessing scripts accepts 3 input parametes: a data file for training (prepared in step 1c), a data file for test (prepared in step 1c) and the directory for output files. Meanwhile, two options are required for specifying the data type (-d) and the directory of structure files (prepared in step 1a) (-s). 
```bash
RME-Preprocess -d SHAPE -s example/dat/structure example/dat/data/SHAPE.train.23SrRNA.data example/dat/data/SHAPE.test.data example/processing-data
RME-Preprocess -d PARS -s example/dat/structure example/dat/data/PARS.train.25SrRNA.data example/dat/data/PARS.test.data example/processing-data
RME-Preprocess -d DMSseq -s example/dat/structure example/dat/data/DMSseq.train.25SrRNA.data example/dat/data/DMSseq.test.data example/processing-data
```

#### b. Output files of step 2

Outfiles contains pairing probabilities and structure information of all bases in the training RNAs (useful for optimizing RME parameters, step 4)
```
example/processing-data/SHAPE.for-optimize-parameter.txt
example/processing-data/PARS.for-optimize-parameter.txt
example/processing-data/DMSseq.for-optimize-parameter.txt
```
Outfiles contains pairing probabilities of all bases in the test RNAs (used for structure prediction, step 3)
```
example/processing-data/SHAPE.for-test.txt
example/processing-data/PARS.for-test.txt
example/processing-data/DMSseq.for-test.txt
```

#### c. Advanced options for RME-Preprocess

We provide several methods for calculating the prior probability of the paired bases, including using a constant (0.535, the default setting), calculating from the training RNAs (train) and calculating from partition function for each RNA (partition). You can select one of them by specify the -r option. Note that when using the partition mode, you should also specify the directory of the sequence files (prepared in step 1b) by -q.

```bash
RME-Preprocess -r constant -d SHAPE -s example/dat/structure example/dat/data/SHAPE.train.23SrRNA.data example/dat/data/SHAPE.test.data example/processing-data
RME-Preprocess -r train -d PARS -s example/dat/structure example/dat/data/PARS.train.25SrRNA.data example/dat/data/PARS.test.data example/processing-data
RME-Preprocess -r partition -q example/dat/sequence -d DMSseq -s example/dat/structure example/dat/data/DMSseq.train.25SrRNA.data example/dat/data/DMSseq.test.data example/processing-data
```

Meanwhile, for SHAPE data, we provided three types of normalization methods, quantile normalization (Bolstad et al., 2003) (quantile, the default setting), 2%/8% rule normalization (Ouyang et al., 2013) (28rule) and no normalization (no). You can select one of them for SHAPE data by specify the -n option.

```bash
RME-Preprocess -n quantile -d SHAPE -s example/dat/structure example/dat/data/SHAPE.train.23SrRNA.data example/dat/data/SHAPE.test.data example/processing-data
RME-Preprocess -n 28rule -d SHAPE -s example/dat/structure example/dat/data/SHAPE.train.23SrRNA.data example/dat/data/SHAPE.test.data example/processing-data
RME-Preprocess -n no -d SHAPE -s example/dat/structure example/dat/data/SHAPE.train.23SrRNA.data example/dat/data/SHAPE.test.data example/processing-data
```
You can get more information about the usage of `RME-Preprocess` by typing `RME-Preprocess -h`.


### Step 3: Structure prediction by RME

#### a. Structure prediction by RME
The RME software accepts 2 input parameters, a data file containing pairing probability for each test RNA (prepared in Step 2b) and a directory for output files. Meanwhile, we recommended to use the -d option to specify the data type, then RME would use default parameter values for that data type. And a -p option specifying the cpu number is quite useful for running in parallel.  

```bash
mkdir example/prediction/
cd example/prediction/
RME -d SHAPE -p 8 ../processing-data/SHAPE.for-test.txt SHAPE 
RME -d PARS -p 8 ../processing-data/PARS.for-test.txt PARS
RME -d DMSseq -p 8 ../processing-data/DMSseq.for-test.txt DMSseq
cd ../..
```

#### b. Output files of the step 3

Outfiles of the predicted structures in CT format
```
example/prediction/SHAPE/*.ct
example/prediction/PARS/*.ct
example/prediction/DMSseq/*.ct
```

#### c. Advanced options for paramter setting

RME allows users to set the parameters in other ways. For example, parameters can be recorded in a file (example: `example/optimize-param/*.params.txt`) (see step 4b for instructions) and fed to RME by using the `-i` option. 
```
mkdir example/prediction/
cd example/prediction/
RME -i ../optimize-param/SHAPE.params.txt -p 8 ../processing-data/SHAPE.for-test.txt SHAPE
RME -i ../optimize-param/PARS.params.txt -p 8 ../processing-data/PARS.for-test.txt PARS
RME -i ../optimize-param/DMSseq.params.txt -p 8 ../processing-data/DMSseq.for-test.txt DMSseq
cd ../..
```

Moreover, the parameters can be set directly by using the `-m`, `--gamma1` and `--gamma2` options. Users are free to adjust the relative contribution of structure probing data by adjusting this values.


#### d. Advanced options for prediction mode

RME combined two rounds of BPPM updating for improved structure prediction, called RME-pre and RME-post, respectively. We provided the `-pre` and `-post` options to enable the users to select either part for structure prediction. This allows the users to test the importance of integrating `RME-pre` and `RME-post` into `RME`.

Moreover, RME-post introduced an `w(i,j)` term to penalize spurious pairs that usually involving in shorter helices. We also provided the `--no-helix-weight` option to enable the users to run RME without this term. This allows the users to test the importance of the term for `RME-post` or `RME`.

You can get more information about usage of `RME` by typing `RME -h`.


### Step 4: optimize the paramter for RME

#### a. Parameter optimization by RME-Optimize

This step is optional. The above 3 sections are sufficient for using RME. This section is for advanced users who want to optimize the RME parameters for their own training data. In this case, try `RME-Optimize`!

`RME-Optimize` accepts two parameters: one is the data file containg pairing probability for the training RNAs (prepared in step 2b), and the other is the output file for recording the optimal parameters for your own data. You can use the `--details` option to specify a file for recording the whole optimizing process. For optimizing in parallel, you can specify the `-p` option to define the cpu number for running. Besides, you can adjust the grid for parameter search by using the `--m-range`, `--gamma1-range` and `--gamma2-range` (providing the min, max and step for grid search seperated by a colon).  

```bash
mkdir example/optimize-param/
cd example/optimize-param/
RME-Optimize ../processing-data/SHAPE.for-optimize-parameter.txt SHAPE.params.txt --m-range 0:1:0.05 --details SHAPE.details.txt -p 4
RME-Optimize ../processing-data/PARS.for-optimize-parameter.txt PARS.params.txt --m-range 0:1:0.02 --details PARS.details.txt -p 4
RME-Optimize ../processing-data/DMSseq.for-optimize-parameter.txt DMSseq.params.txt --m-range 0:1:0.02 --details DMSseq.details.txt -p 4
cd ../..
```

You can get more information about usage of `RME-Optimize` by typing `RME-Optimize -h`.

#### b. Output files of step 4

Output files containing the optimal parameters, which can be used in step 3c
```
example/optimize-param/SHAPE.params.txt
example/optimize-param/PARS.params.txt
example/optimize-param/DMSseq.params.txt
```

Output files recording the detail of the optimizing processes
```
example/optimize-param/SHAPE.details.txt
example/optimize-param/PARS.details.txt
example/optimize-param/DMSseq.details.txt
```

## Output files

structure prediction files

```
example/prediction/*/*.ct
```

pairing probability files

```
example/processing-data/*/*.prob
```

