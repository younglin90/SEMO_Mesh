CCOMPLR = mpiicpc
CFLAGS = -c -O3
EXE = SEMO
LIBINCLUDE = -I./lib/Metis\
             -I./lib/Scotch\

LIBRARIES = ./lib/Metis/libparmetis.a\
            ./lib/Metis/libmetis.a\
            ./lib/Scotch/libscotcherr.a\
            ./lib/Scotch/libscotch.a\

SOURCES = src/main.cpp\
	      src/mesh/build.cpp\
	      src/mesh/load.cpp\
	      src/mesh/save.cpp\
	      src/mesh/geometric.cpp\
	      src/mesh/partition.cpp\
	      src/mesh/distribute.cpp\
	      src/math/functions.cpp\
	      src/mpi/build.cpp\

OBJECTS = $(SOURCES:.cpp=.o)

all : $(EXE)

$(EXE) : $(OBJECTS)
	$(CCOMPLR) -o $@ $(OBJECTS) $(LIBRARIES)

%.o : %.cpp
	$(CCOMPLR) $(CFLAGS) $(LIBINCLUDE) $< -o $@

clean:
	@echo -e "\033[1;31m deleting objects \033[0m" | tee make.log
	@rm -fr $(OBJECTS) $(EXE) *.o