# complier options


PREFIX = /usr/local
ROOTPATH = .
RSPATH = RNAstructure-src
BINPATH = $(PREFIX)/bin
EXEPATH = $(PREFIX)/libexec/RME
DATPATH = $(PREFIX)/share/RME/data_tables

CXX = g++
CXXFLAGS = -O3 -fsched-spec-load -I${RSPATH}
LDFLAGS = -s
LINK = ${CXX} ${LDFLAGS} -o $@
INSTALL = install

OUTDIRS = bin bin/exec

PROC_EXE = distributionFitting.R getFullLenAndStruct.pl \
	posteriorPairingProb.R quantileNormalization.R \
	readCountNormalization.R
PROC_EXE_IN = $(addprefix processing/, $(PROC_EXE))
PROC_EXE_OUT = $(addprefix bin/exec/, $(PROC_EXE))
PROC_EXE_INST = $(addprefix $(EXEPATH)/, $(PROC_EXE))

PROC_BIN = processDMSseq processSHAPE processPARS
PROC_BIN_IN = $(addsuffix .in, $(addprefix processing/, $(PROC_BIN)))
PROC_BIN_OUT = $(addprefix bin/, $(PROC_BIN))

ALL_BIN = RME RME-Optimize processSHAPE processDMSseq processPARS

# The text interface command line parser.
CMD_LINE_PARSER = \
	${RSPATH}/src/ParseCommandLine.o
# The utility that scores a structure against another.
STRUCTURE_SCORER = \
	${RSPATH}/src/score.o
# Common files for the RNA library.
RNA_FILES = \
	${RSPATH}/RNA_class/RNA.o \
	${RSPATH}/RNA_class/thermodynamics.o \
	${RSPATH}/src/algorithm.o \
	${RSPATH}/src/alltrace.o \
	${RSPATH}/src/arrayclass.o \
	${RSPATH}/src/dotarray.o \
	${RSPATH}/src/draw.o \
	${RSPATH}/src/extended_double.o \
	${RSPATH}/src/forceclass.o \
	${RSPATH}/src/MaxExpect.o \
	${RSPATH}/src/MaxExpectStack.o \
	${RSPATH}/src/outputconstraints.o \
	${RSPATH}/src/pfunction.o \
	${RSPATH}/src/probknot.o \
	${RSPATH}/src/random.o \
	${RSPATH}/src/rna_library.o \
	${RSPATH}/src/stackclass.o \
	${RSPATH}/src/stackstruct.o \
	${RSPATH}/src/stochastic.o \
	${RSPATH}/src/structure.o \
	${RSPATH}/src/TProgressDialog.o

RME_FILES = \
	${ROOTPATH}/RME/RMEPre.o \
	${ROOTPATH}/RME/RMEPost.o \
	${ROOTPATH}/RME/Utils.o

.PHONY: RME scorer ct2dot dot2ct processing

all: $(OUTDIRS) RME processing scorer ct2dot dot2ct

install: all
	# create target directory
	$(INSTALL) -d ${BINPATH} $(EXEPATH) ${DATPATH}
	# install programs
	$(INSTALL) -T RME/RME.in ${BINPATH}/RME
	$(INSTALL) -T RME/RME-Optimize.in $(BINPATH)/RME-Optimize
	$(INSTALL) -T processing/processSHAPE.in $(BINPATH)/processSHAPE
	$(INSTALL) -T processing/processPARS.in $(BINPATH)/processPARS
	$(INSTALL) -T processing/processDMSseq.in $(BINPATH)/processDMSseq
	# substitute variables
	sed -i "s#{DATAPATH}#'$(DATPATH)'#;s#{EXEPATH}#'$(EXEPATH)'#" $(BINPATH)/RME
	sed -i "s#{DATAPATH}#'$(DATPATH)'#;s#{EXEPATH}#'$(EXEPATH)'#" $(BINPATH)/RME-Optimize
	sed -i "s#{EXEPATH}#'$(EXEPATH)'#" $(BINPATH)/processSHAPE
	sed -i "s#{EXEPATH}#'$(EXEPATH)'#" $(BINPATH)/processPARS
	sed -i "s#{EXEPATH}#'$(EXEPATH)'#" $(BINPATH)/processDMSseq
	# install required data files
	$(INSTALL) -t $(EXEPATH) bin/exec/*
	$(INSTALL) --mode 644 -t $(DATPATH) data_tables/*
	
subvars-inst: $(addprefix $(PREFIX)/$(BINPATH), $(ALL_BIN))
	sed -i "s#{DATAPATH}#'$(DATPATH)'#; s#{EXEPATH}#'$(EXEPATH)'#" $^
	
$(OUTDIRS):
	mkdir -p $@
	
# Build the ct2dot text interface.
ct2dot: bin/ct2dot
bin/ct2dot: ${RSPATH}/ct2dot/ct2dot.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} ${RSPATH}/ct2dot/ct2dot.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the dot2ct text interface.
dot2ct: bin/dot2ct
bin/dot2ct: ${RSPATH}/dot2ct/dot2ct.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} ${RSPATH}/dot2ct/dot2ct.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the scorer interface.
scorer: bin/scorer
bin/scorer: ${RSPATH}/scorer/Scorer_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}
	${LINK} ${RSPATH}/scorer/Scorer_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}
	
# Build RME
RME: bin/RME bin/RME-Optimize
	
bin/RME: RME/RME.in bin/exec/RME 
	sed "s#{DATAPATH}#'${CURDIR}/data_tables'#; s#{EXEPATH}#'${CURDIR}/bin/exec'#" $< > $@
	chmod 755 $@
	
bin/RME-Optimize:  RME/RME-Optimize.in bin/exec/RME-Optimize
	sed "s#{DATAPATH}#'${CURDIR}/data_tables'#; s#{EXEPATH}#'${CURDIR}/bin/exec'#" $< > $@
	chmod 755 $@
	
bin/exec/RME:  RME/RMEInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} 
	${LINK} ${ROOTPATH}/RME/RMEInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} -lpthread
	
bin/exec/RME-Optimize:  RME/RMEOptimizeInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} 
	${LINK} ${ROOTPATH}/RME/RMEOptimizeInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} -lpthread 

processing: $(PROC_BIN_OUT) $(PROC_EXE_OUT)

bin/processSHAPE: processing/processSHAPE.in
	sed "s#{EXEPATH}#'${CURDIR}/bin/exec'#" $< > $@
	chmod 755 $@
	
bin/processPARS: processing/processPARS.in
	sed "s#{EXEPATH}#'${CURDIR}/bin/exec'#" $< > $@
	chmod 755 $@
	
bin/processDMSseq: processing/processDMSseq.in
	sed "s#{EXEPATH}#'${CURDIR}/bin/exec'#" $< > $@
	chmod 755 $@
	
$(PROC_EXE_OUT):
	cp $(PROC_EXE_IN) bin/exec
	
%.o: %.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<
##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables
##########

# Remove object files and any temporary files from building.
clean:
	find . -depth -name '*~' -delete
	find . -depth -name '*.o' -delete
	rm -rf bin

