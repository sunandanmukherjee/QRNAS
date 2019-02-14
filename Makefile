CXX = g++
CXXFLAGS = --std=c++03 -Wall -O3 -Iinclude
LDFLAGS_P = -pthread
LDFLAGS_S = -DSEQUENTIAL

SOURCES = bprestr.cpp main.cpp tatom.cpp tbond.cpp tdihed.cpp textutils.cpp timpro.cpp tmolecule.cpp tposrestr.cpp tvector.cpp distrestr.cpp tangle.cpp tbasepair.cpp tconsts.cpp telectr.cpp thbond.cpp tlennard.cpp tnonbond.cpp tspring.cpp

OBJECTS = $(SOURCES:.cpp=.o)


TARGET = QRNA

default: sequential

sequential:
	$(CXX) $(CXXFLAGS) $(LDFLAGS_S) -o $(TARGET) $(SOURCES)
#	$(CXX) $(CXXFLAGS) $(LDFLAGS_S) -c $(SOURCES)
#	$(CXX) $(CXXFLAGS) $(LDFLAGS_S) -o $(TARGET) $(OBJECTS)

parallel:
	$(CXX) $(CXXFLAGS) $(LDFLAGS_P) -o $(TARGET) $(SOURCES)
#	$(CXX) $(CXXFLAGS) $(LDFLAGS_P) -c $(SOURCES)
#	$(CXX) $(CXXFLAGS) $(LDFLAGS_P) -o $(TARGET) $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(TARGET)
# 
