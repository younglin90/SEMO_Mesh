CCOMPLR = mpiicpc
CFLAGS = -c -O3
EXE = SEMO
SOURCES = src/main.cpp\
	      src/mesh/meshLoad.cpp\
	      src/mesh/meshConnector.cpp\
	      src/mesh/meshGeometric.cpp\

OBJECTS = $(SOURCES:.cpp=.o)

all : $(EXE)

$(EXE) : $(OBJECTS)
	$(CCOMPLR) -o $@ $(OBJECTS)

%.o : %.cpp
	$(CCOMPLR) $(CFLAGS) $< -o $@
