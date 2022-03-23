#Gnu makefile for gnu+Linux system
CC=g++

IDIR =include
ODIR=obj
LOGDIR=log
SRCDIR=src
OUTDIR=bin
DATADIR=out
OUTNAME=example.exe

CXXFLAGS=-std=c++17 -Wall -O2 -Wextra -Wpedantic -Wdouble-promotion -I$(IDIR)

LIBS = -lSDL2 -lSDL2_image -lSDL2_mixer -lGLEW -lGL

_OBJ = main.o field.o constants.o particle.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))



#I can't place literal tabs, tab default to space x4. This seems better
.RECIPEPREFIX := ;


#Compile the final program
$(OUTDIR)/$(OUTNAME):	$(OBJ)
;$(CC) -o $@ $^ $(CXXFLAGS)


#Create object files, I write these individually, since I don't know a way to automate dependency detection, the alternative would be to make all .hpp files dependencies of everything, but that would be impractical
$(ODIR)/main.o: $(SRCDIR)/main.cpp $(IDIR)/field.hpp  $(IDIR)/constants.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/field.o: $(SRCDIR)/field.cpp $(IDIR)/field.hpp  $(IDIR)/constants.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/constants.o: $(SRCDIR)/constants.cpp $(IDIR)/constants.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/particle.o: $(SRCDIR)/particle.cpp $(IDIR)/constants.hpp $(IDIR)/field.hpp  $(IDIR)/particle.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS)

$(OUTDIR)/test_double_precision.exe: $(SRCDIR)/test_double_precision.cpp
;$(CC) -o $@ $< $(CXXFLAGS)

.PHONY: clean
.PHONY: clean_log
.PHONY: clean_out

clean_log:
;rm -f $(LOGDIR)/*


clean_out:
;rm -fr $(DATADIR)/out*

clean:
;rm -f $(ODIR)/*.o
;rm -f $(OUTDIR)/$(OUTNAME)

