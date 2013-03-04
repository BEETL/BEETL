CFLAGS:=-O3 -I. -Wall
SEQAN_DIR:=/home/acox/Indel/tony/BEETL/Seqan1.3/
#SEQAN_DIR:=/illumina/development/SeqAn/
TARGETS:=BCRext

#Timer.o: Timer.cpp
#	g++ -c $(CFLAGS) Timer.cpp

#BCRext.o: BCRext.cpp BCRext.hh
#	 g++ -c $(CFLAGS) BCRext.cpp
#ReadBuffer.o: ReadBuffer.cpp ReadBuffer.hh
#	 g++ -c $(CFLAGS) ReadBuffer.cpp

%.o: %.cpp %.hh
	g++ -c $(CFLAGS) $<

BCRext: Timer.o ReadBuffer.o BwtReader.o BCRext.o
	g++ $(CFLAGS) -o $@ $^

countWords: Timer.o ReadBuffer.o BwtReader.o CountWords.o
	g++ $(CFLAGS) -o $@ $^
	
makeBWT: Timer.o MakeBWT.o
	g++ $(CFLAGS) -o $@ $^

makeBWTSkew: BuildChromosomeBwt.cpp
	c++ $(CFLAGS) -I $(SEQAN_DIR) -o $@ $^

mergeBWT: MergeChromosomeBwt.cpp
	g++ $(CFLAGS) -o $@ $^


all: ${TARGETS}

clean: /bin/rm -f *.o *~ ${TARGETS}

release:
	mkdir BEETL
	cp Alphabet.hh BwtReader.hh ReadBuffer.hh Types.hh BCRext.cpp BwtWriter.hh LetterCount.hh Timer.cpp BCRext.hh Config.hh Timer.hh BwtReader.cpp ReadBuffer.cpp makefile BEETL
	cp Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf BEETL
	tar -cjvf BEETL.tar.bz2 BEETL/
	@echo release tar file made successfully 
