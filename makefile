CC=g++
CXXFILES = YRoots/main.cpp
CXXTESTFILES = TestYRoots/TestMain.mm

IDIR = -IYRoots/include -IYRoots/3rdParty -IYRoots/FFTW/include -I/usr/local/include
LDIR = -L/usr/local/include -LYRoots/FFTW/lib
IDIRTEST = -ITestYRoots

LIBS = -lfftw3 -lm -pthread
CXXFLAGS = -std=c++11 -Wall -O3
MACROS = -Wextra -DUSE_TIMING

all:
	$(CC) -o yroots_solver $(CXXFLAGS) $(IDIR) $(LDIR) $(CXXFILES) $(LIBS) $(MACROS)

test:
	$(CC) -o yroots_test $(CXXFLAGS) $(IDIR) $(IDIRTEST) $(LDIR) $(CXXTESTFILES) $(LIBS) $(MACROS)