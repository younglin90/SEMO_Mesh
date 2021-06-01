CCOMPLR = mpiicpc

CFLAGS = -c -O3

EXE = SEMO

LIBINCLUDE = \
             -I./lib/zlib\
			 -I/home/yyl/petsc/include\
             -I./lib/PETSc/include\
             -I./lib/Metis\
             -I./lib/Scotch\
			 # -I/home/yyl/petsc/arch-linux-c-debug/include\

LIBRARIES = \
            ./lib/zlib/libz.a\
			./lib/PETSc/lib/libpetsc.so.3.14\
			-Wl,-rpath,./lib/PETSc/lib\
            ./lib/Metis/libparmetis.a\
            ./lib/Metis/libmetis.a\
            ./lib/Scotch/libscotch.a\
            ./lib/Scotch/libscotcherr.a\
            # ./lib/Scotch/libscotchmetis.a\
			# -Wl,-rpath,/home/yyl/petsc/arch-linux-c-debug/lib\
			# -L/home/yyl/petsc/arch-linux-c-debug/lib\
			# -lm\
			# -lpetsc\
			# ./lib/PETSc/lib/libpetsc.so.3.14\
			# -Wl,-rpath,./lib/PETSc/lib\

SOURCES = src/main.cpp\
	      src/mesh/build.cpp\
	      src/mesh/load.cpp\
	      src/mesh/save.cpp\
	      src/mesh/geometric.cpp\
	      src/mesh/partition.cpp\
	      src/mesh/distribute.cpp\
	      src/mesh/hexaOctAMR.cpp\
	      src/math/math.cpp\
	      src/math/gradient.cpp\
	      src/mpi/build.cpp\
	      src/controls/build.cpp\
	      src/variables/build.cpp\
	      src/solvers/build.cpp\
	      src/solvers/pressureBased.cpp\
	      src/solvers/densityBased.cpp\
	      src/solvers/reconIncom.cpp\
	      src/solvers/reconComp.cpp\
	      src/solvers/eqMomentum.cpp\
	      src/solvers/eqPressure.cpp\
	      src/solvers/eqVolfrac.cpp\
	      src/solvers/solvePETSc.cpp\
	      src/solvers/timestep.cpp\
	      src/solvers/RHS.cpp\
	      src/solvers/flux.cpp\
	      src/solvers/source.cpp\
	      src/solvers/linearSolver.cpp\
	      src/solvers/eos.cpp\
	      src/solvers/norm.cpp\
	      src/solvers/NVD.cpp\
	      src/solvers/hybridBased.cpp\
	      src/utility/read.cpp\

OBJECTS = $(SOURCES:.cpp=.o)

EXE_PARTITION = PARTITION

SOURCES_PARTITION = src/utility/partition.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/utility.cpp\

OBJECTS_PARTITION = $(SOURCES_PARTITION:.cpp=.o)

COTEXT  = "\033[1;31m Compiling\033[0m\033[1m $< \033[0m"

EXE_LOAD = LOAD

SOURCES_LOAD = src/utility/load.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/utility.cpp\

OBJECTS_LOAD = $(SOURCES_LOAD:.cpp=.o)

EXE_INITIAL = INITIAL

SOURCES_INITIAL = src/utility/initial.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/utility.cpp\
                    src/controls/build.cpp\
                    src/utility/read.cpp\

OBJECTS_INITIAL = $(SOURCES_INITIAL:.cpp=.o)

COTEXT  = "\033[1;31m Compiling\033[0m\033[1m $< \033[0m"

# all : $(EXE)
all : $(EXE) $(EXE_PARTITION) $(EXE_INITIAL)
# all : $(EXE_LOAD)

$(EXE) : $(OBJECTS)
	@$(CCOMPLR) -o $@ $(OBJECTS) $(LIBRARIES)
	@echo " "
	@echo "\033[1;31m CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_PARTITION) : $(OBJECTS_PARTITION)
	@$(CCOMPLR) -o $@ $(OBJECTS_PARTITION) $(LIBRARIES)
	@echo " "
	@echo "\033[1;31m CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_LOAD) : $(OBJECTS_LOAD)
	@$(CCOMPLR) -o $@ $(OBJECTS_LOAD) $(LIBRARIES)
	@echo " "
	@echo "\033[1;31m CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_INITIAL) : $(OBJECTS_INITIAL)
	@$(CCOMPLR) -o $@ $(OBJECTS_INITIAL) $(LIBRARIES)
	@echo " "
	@echo "\033[1;31m CODE compile/link complete \033[0m" | tee -a make.log

%.o : %.cpp
	@echo $(COTEXT) | tee -a make.log
	@$(CCOMPLR) $(CFLAGS) $(LIBINCLUDE) $< -o $@

clean:
	@echo "\033[1;31m deleting objects \033[0m" | tee make.log
	@rm -fr $(OBJECTS) $(EXE) $(EXE_PARTITION) $(EXE_INITIAL) make.log *.o