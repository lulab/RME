# complier options


PREFIX = /usr/local
ROOTPATH = .
RSPATH = RNAstructure-src
BINPATH = bin
EXEPATH = libexec/RME
DATAPATH = share/RME/data_tables

CXX = g++
CXXFLAGS = -O2 -I${RSPATH}
LDFLAGS = -s
LINK = ${CXX} ${LDFLAGS} -o $@
INSTALL = install

OUTDIRS = bin libexec/RME

PROC_EXE = distributionFitting.R getFullLenAndStruct.pl \
	posteriorPairingProb.R quantileNormalization.R \
	readCountNormalization.R 28ruleNormalization.R \
	processDMSseq processSHAPE processPARS
PROC_EXE_IN = $(addprefix processing/, $(PROC_EXE))
PROC_EXE_OUT = $(addprefix libexec/RME/, $(PROC_EXE))
PROC_EXE_INST = $(addprefix $(EXEPATH)/, $(PROC_EXE))

PROC_BIN = RME-Preprocess
PROC_BIN_IN = $(addsuffix .in, $(addprefix processing/, $(PROC_BIN)))
PROC_BIN_OUT = $(addprefix bin/, $(PROC_BIN))

ALL_BIN = RME RME-Optimize RME-Preprocess

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

.PHONY: RME processing

all: $(OUTDIRS) data_tables RME processing scorer ct2dot dot2ct partition 

install: all
	# create target directory
	$(INSTALL) -d $(PREFIX)/${BINPATH} $(PREFIX)/$(EXEPATH) $(PREFIX)/${DATAPATH}
	# install programs
	$(INSTALL) -T RME/RME.in $(PREFIX)/${BINPATH}/RME
	$(INSTALL) -T RME/RME-Optimize.in $(PREFIX)/$(BINPATH)/RME-Optimize
	$(INSTALL) -T processing/RME-Preprocess.in $(PREFIX)/$(BINPATH)/RME-Preprocess
	# substitute variables
	sed -i "s#{DATAPATH}#$(DATAPATH)#;s#{EXEPATH}#$(EXEPATH)#" $(PREFIX)/$(BINPATH)/RME
	sed -i "s#{DATAPATH}#$(DATAPATH)#;s#{EXEPATH}#$(EXEPATH)#" $(PREFIX)/$(BINPATH)/RME-Optimize
	sed -i "s#{DATAPATH}#$(DATAPATH)#;s#{EXEPATH}#$(EXEPATH)#" $(PREFIX)/$(BINPATH)/RME-Preprocess
	# install required data files
	$(INSTALL) -t $(PREFIX)/$(EXEPATH) $(EXEPATH)/*
	$(INSTALL) --mode 644 -t $(PREFIX)/$(DATAPATH) $(DATAPATH)/*
	
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
	
partition: bin/partition
${EXEPATH}/partition: ${RSPATH}/pfunction/partition.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}
	${LINK} ${RSPATH}/pfunction/partition.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}
	
# Build RME
RME: bin/RME bin/RME-Optimize bin/partition ${EXEPATH}/AveragePairProb
	
bin/RME: RME/RME.in ${EXEPATH}/RME 
	sed "s#{DATAPATH}#${DATAPATH}#; s#{EXEPATH}#${EXEPATH}#" $< > $@
	chmod 755 $@
	
bin/RME-Optimize:  ${ROOTPATH}/RME/RME-Optimize.in ${EXEPATH}/RME-Optimize
	sed "s#{DATAPATH}#${DATAPATH}#; s#{EXEPATH}#${EXEPATH}#" $< > $@
	chmod 755 $@
	
bin/partition:  ${ROOTPATH}/RME/partition.in ${EXEPATH}/partition
	sed "s#{DATAPATH}#${DATAPATH}#; s#{EXEPATH}#${EXEPATH}#" $< > $@
	chmod 755 $@
	
${EXEPATH}/RME:  ${ROOTPATH}/RME/RMEInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} 
	${LINK} ${ROOTPATH}/RME/RMEInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} -lpthread
	
${EXEPATH}/RME-Optimize:  ${ROOTPATH}/RME/RMEOptimizeInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} 
	${LINK} ${ROOTPATH}/RME/RMEOptimizeInterface.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} -lpthread 

${EXEPATH}/AveragePairProb:  ${ROOTPATH}/RME/AveragePairProb.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} 
	${LINK} ${ROOTPATH}/RME/AveragePairProb.o ${RME_FILES} ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES} -lpthread

processing: $(PROC_EXE_OUT) $(PROC_BIN_OUT)
	
data_tables: $(DATAPATH)

$(DATAPATH): $(RSPATH)/data_tables
	mkdir -p $(DATAPATH)
	cp $</* $@
	
bin/RME-Preprocess: $(ROOTPATH)/processing/RME-Preprocess.in
	sed "s#{DATAPATH}#${DATAPATH}#; s#{EXEPATH}#${EXEPATH}#" $< > $@
	chmod 755 $@
	
$(PROC_EXE_OUT):
	cp $(PROC_EXE_IN) ${EXEPATH}
	
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
	rm -rf bin share libexec

