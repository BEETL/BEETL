#
# Copyright (c) 2011 Illumina, Inc.
#
# 
# This software is covered by the "Illumina Non-Commercial Use Software
# and Source Code License Agreement" and any user of this software or
# source file is bound by the terms therein (see accompanying file
# Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
#
# This file is part of the BEETL software package.
#
# Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
# Lightweight BWT Construction for Very Large String Collections. 
# Proceedings of CPM 2011, pp.219-231
#

# Makefile for the BEETL software library
#
# Makefile version $Id: Makefile,v 1.7 2011/11/28 16:38:32 acox Exp $


# >> SETUP SECTION

# compiler: e.g. g++
CC=g++

# include path for header files
INC:= -I. -Ishared/ -IBCR/ -IBCRext/  

# DEBUG build build with make DEBUG
ifdef DEBUG
CFLAGS:= $(INC) -Wall -g #-D DEBUG
else
# normal build, optimized for speed
CFLAGS:= $(INC) -Wall -O3
endif

TARGETS:= Beetl

BEETL_DIR_NAME := beetl/	

OBJ_DIR:= ../build

####### normally you don't need to change anything below this line ########

### main target, Beetl executable	

Beetl: Beetl.o BCRext.o BCR_BWTCollection.o Timer.o Tools.o\
       BWTCollection.o BCRexternalBWT.o BuildBCR.o Sorting.o TransposeFasta.o BwtReader.o BwtWriter.o \
	ReadBuffer.o 
	cd ../build/; ${CC} $(CFLAGS) -o ../$@ $^

### shared stuff

Beetl.o: shared/Beetl.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@ 

Tools.o: shared/Tools.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@

Timer.o: shared/Timer.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@

### BCRext

ReadBuffer.o: BCRext/ReadBuffer.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@	

BCRext.o: BCRext/BCRext.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@	

BwtReader.o: BCRext/BwtReader.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@

BwtWriter.o: BCRext/BwtWriter.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@

### BCR	

BWTCollection.o: BCR/BWTCollection.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@		

BCR_BWTCollection.o: BCR/BCR_BWTCollection.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@		

TransposeFasta.o: BCR/TransposeFasta.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@	

Sorting.o: BCR/Sorting.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@			

Transpose.o: BCR/Transpose.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@	

BCRexternalBWT.o: BCR/BCRexternalBWT.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@		

BuildBCR.o: BCR/BuildBCR.cpp
	${CC} -c $(CFLAGS) $< -o $(OBJ_DIR)/$@		



### >> done compiling / linking

## special targets

all: ${TARGETS}
	@echo BEETL compiled successfully

# normal clean, build files and binary
clean: 
	@echo cleaning up build and executable files
	/bin/rm -f *~ ${TARGETS}
	/bin/rm -f ../build/*.o
	/bin/rm -f ../Beetl	

# cleans also files generated during "make release"
cleaner: clean
	@echo cleaning up release files
	cd ../..; /bin/rm -f BEETL.tar.bz2 

###### release target, no files in script/ + all CVS files are skipped

release: all cleaner	    
	    cd ../..; tar -cjvf BEETL.tar.bz2 ${BEETL_DIR_NAME} --exclude=*CVS* --exclude=scripts/*
	    @echo release tar file made successfully
