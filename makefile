CC=g++
CXXFILES = YRoots/main.cpp
CXXTESTFILES = TestYRoots/TestMain.mm

IDIR = -IYRoots/include -IYRoots/3rdParty
LDIR = -LYRoots/lib
IDIRTEST = -ITestYRoots

LIBS = -lfftw3 -lm
CXXFLAGS = -std=c++11 -Wall -Ofast
MACROS = -Wextra -DUSE_TIMING

all:
	$(CC) -o yroots_solver $(CXXFLAGS) $(IDIR) $(LDIR) $(CXXFILES) $(LIBS) $(MACROS)

test:
	$(CC) -o yroots_test $(CXXFLAGS) $(IDIR) &(IDIRTEST) $(LDIR) $(CXXTESTFILES) $(LIBS) $(MACROS)