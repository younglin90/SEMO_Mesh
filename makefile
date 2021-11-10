CCOMPLR = mpiicpc

CFLAGS = -std=c++17 -c -O2

EXE = SEMO

LIBINCLUDE = \
             -Ilib/zlib\
             -Ilib/PETSc/include\
             -Ilib/Metis\
             -Ilib/Scotch\
             -Ilib/HYPRE/include\
             -Ilib/amgcl\
             -Ilib/boost\
             # -I/home/yyl/petsc/arch-linux-c-debug/include\

LIBRARIES = \
            lib/zlib/libz.a\
            lib/PETSc/lib/libpetsc.so.3.14\
            -Wl,-rpath,lib/PETSc/lib\
            lib/Metis/libparmetis.a\
            lib/Metis/libmetis.a\
            lib/Scotch/libscotch.a\
            lib/Scotch/libscotcherr.a\
            lib/HYPRE/lib/libHYPRE.a\
             # ./lib/Scotch/libscotchmetis.a\
             # -Wl,-rpath,/home/yyl/petsc/arch-linux-c-debug/lib\
             # -L/home/yyl/petsc/arch-linux-c-debug/lib\
             # -lm\
             # -lpetsc\
             # ./lib/PETSc/lib/libpetsc.so.3.14\
             # -Wl,-rpath,./lib/PETSc/lib\

SOURCES = src/mesh/build.cpp\
	      src/mesh/load.cpp\
	      src/mesh/save.cpp\
	      src/mesh/geometric.cpp\
	      src/mesh/partition.cpp\
	      src/mesh/distribute.cpp\
	      src/mesh/hexaOctAMR.cpp\
	      src/math/math.cpp\
	      src/math/gradient.cpp\
	      src/solvers/transport.cpp\
	      src/solvers/viscousFlux.cpp\
	      src/mesh/saveGnuplot.cpp\
	      src/mpi/build.cpp\
	      src/controls/build.cpp\
	      src/variables/build.cpp\
	      src/solvers/build.cpp\
	      src/solvers/pressureBased.cpp\
	      src/solvers/densityBased.cpp\
	      src/solvers/compCoupled.cpp\
	      src/solvers/compressible/massfrac.cpp\
	      src/solvers/compressible/pressure.cpp\
	      src/solvers/compressible/momentum.cpp\
	      src/solvers/compressible/energy.cpp\
	      src/solvers/compressible/flows.cpp\
	      src/solvers/compressible/coupled.cpp\
	      src/solvers/reconIncom.cpp\
	      src/solvers/reconComp.cpp\
	      src/solvers/eqCoupled.cpp\
	      src/solvers/eqMomentum.cpp\
	      src/solvers/eqPressure.cpp\
	      src/solvers/eqVolfrac.cpp\
	      src/solvers/solveAMGCL.cpp\
	      src/solvers/solvePETSc.cpp\
	      src/solvers/solveHYPRE.cpp\
	      src/solvers/timestep.cpp\
	      src/solvers/RHS.cpp\
	      src/solvers/flux.cpp\
	      src/solvers/source.cpp\
	      src/solvers/linearSolver.cpp\
	      src/solvers/eos.cpp\
	      src/solvers/norm.cpp\
	      src/solvers/NVD.cpp\
	      src/solvers/hybridBased.cpp\
	      src/solvers/curvature.cpp\
	      src/utility/read.cpp\
	      src/turbulenceModels/LES.cpp\
	      src/math/RCM.cpp\
	      src/mesh/polyAMR.cpp\
	      src/mesh/polyRefine.cpp\
	      src/mesh/polyUnrefine.cpp\
	      src/mesh/reorder.cpp\

OBJECTS = src/main.o $(SOURCES:.cpp=.o)

# Density based single time
EXE_CompDensitySingle = CompDensitySingle

OBJECTS_CompDensitySingle = src/main/compressible/densityBasedSingle.o $(SOURCES:.cpp=.o)

# Density based dual time
EXE_CompDensityDual = CompDensityDual

OBJECTS_CompDensityDual = src/main/compressible/densityBasedDual.o $(SOURCES:.cpp=.o)

# Pressure based 
EXE_IncomPressure = IncomPressure

OBJECTS_IncomPressure = src/main/incompressible/pressureBased.o $(SOURCES:.cpp=.o)

# hybrid based
EXE_CompHybrid = CompHybrid

OBJECTS_CompHybrid = src/main/compressible/hybridBased.o $(SOURCES:.cpp=.o)

# comp coupled based
EXE_CompCoupled = CompCoupled

OBJECTS_CompCoupled = src/main/compressible/mainCoupled.o $(SOURCES:.cpp=.o)

# partitioning
EXE_PARTITION = Partition

SOURCES_PARTITION = src/utility/partition.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/math.cpp\

OBJECTS_PARTITION = $(SOURCES_PARTITION:.cpp=.o)

# initialization
EXE_INITIAL = Initial

SOURCES_INITIAL = src/utility/initial.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/math.cpp\
                    src/controls/build.cpp\
                    src/utility/read.cpp\

OBJECTS_INITIAL = $(SOURCES_INITIAL:.cpp=.o)

# potential flow
EXE_POTENTIAL = Potential

SOURCES_POTENTIAL = src/utility/potential.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/math.cpp\
                    src/math/gradient.cpp\
                    src/mpi/build.cpp\
                    src/controls/build.cpp\
                    src/utility/read.cpp\
                    src/solvers/build.cpp\
                    src/solvers/eos.cpp\
                    src/solvers/solvePETSc.cpp\
                    src/solvers/reconIncom.cpp\
                    src/solvers/reconComp.cpp\
                    src/solvers/NVD.cpp\

OBJECTS_POTENTIAL = $(SOURCES_POTENTIAL:.cpp=.o)

# MapField
EXE_MapField = MapField

SOURCES_MapField = src/utility/mapField.cpp\
                    src/mesh/build.cpp\
                    src/mesh/load.cpp\
                    src/mesh/save.cpp\
                    src/mesh/geometric.cpp\
                    src/math/math.cpp\
                    src/controls/build.cpp\
                    src/utility/read.cpp\

OBJECTS_MapField = $(SOURCES_MapField:.cpp=.o)

# MeshAMR
EXE_MeshAMR = MeshAMR

# SOURCES_MeshAMR = src/utility/meshAMR.cpp\
                    # src/mesh/build.cpp\
                    # src/mesh/load.cpp\
                    # src/mesh/save.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/controls/build.cpp\
                    # src/utility/read.cpp\

OBJECTS_MeshAMR = src/utility/meshAMR.o $(SOURCES:.cpp=.o)
# OBJECTS_MeshAMR = $(SOURCES_MeshAMR:.cpp=.o)

# ExtractData
EXE_ExtractData = ExtractData

OBJECTS_ExtractData = src/utility/extractData.o $(SOURCES:.cpp=.o)

COTEXT  = "\033[1;31m Compiling\033[0m\033[1m $< \033[0m"

# all : $(EXE)
all : $(EXE) $(EXE_CompDensitySingle) $(EXE_CompDensityDual) $(EXE_IncomPressure) $(EXE_CompHybrid) $(EXE_CompCoupled) $(EXE_PARTITION) $(EXE_INITIAL) $(EXE_POTENTIAL) $(EXE_MapField) $(EXE_MeshAMR) $(EXE_ExtractData)
# all : $(EXE_LOAD)

$(EXE) : $(OBJECTS)
	@$(CCOMPLR) -o $@ $(OBJECTS) $(LIBRARIES)
	@echo -e "\033[1;31m Main code compile/link complete \033[0m" | tee -a make.log

$(EXE_CompDensitySingle) : $(OBJECTS_CompDensitySingle)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompDensitySingle) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_density_single CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompDensityDual) : $(OBJECTS_CompDensityDual)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompDensityDual) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_density_dual CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_IncomPressure) : $(OBJECTS_IncomPressure)
	@$(CCOMPLR) -o $@ $(OBJECTS_IncomPressure) $(LIBRARIES)
	@echo -e "\033[1;31m Incom_pressure CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompHybrid) : $(OBJECTS_CompHybrid)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompHybrid) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_hybrid CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompCoupled) : $(OBJECTS_CompCoupled)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompCoupled) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_Coupled CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_PARTITION) : $(OBJECTS_PARTITION)
	@$(CCOMPLR) -o $@ $(OBJECTS_PARTITION) $(LIBRARIES)
	@echo -e "\033[1;31m PARTITION CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_INITIAL) : $(OBJECTS_INITIAL)
	@$(CCOMPLR) -o $@ $(OBJECTS_INITIAL) $(LIBRARIES)
	@echo -e "\033[1;31m INITIAL CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_POTENTIAL) : $(OBJECTS_POTENTIAL)
	@$(CCOMPLR) -o $@ $(OBJECTS_POTENTIAL) $(LIBRARIES)
	@echo -e "\033[1;31m POTENTIAL CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_MapField) : $(OBJECTS_MapField)
	@$(CCOMPLR) -o $@ $(OBJECTS_MapField) $(LIBRARIES)
	@echo -e "\033[1;31m MapField CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_MeshAMR) : $(OBJECTS_MeshAMR)
	@$(CCOMPLR) -o $@ $(OBJECTS_MeshAMR) $(LIBRARIES)
	@echo -e "\033[1;31m MeshAMR CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_ExtractData) : $(OBJECTS_ExtractData)
	@$(CCOMPLR) -o $@ $(OBJECTS_ExtractData) $(LIBRARIES)
	@echo -e "\033[1;31m ExtractData CODE compile/link complete \033[0m" | tee -a make.log

%.o : %.cpp
	@echo -e $(COTEXT) | tee -a make.log
	@$(CCOMPLR) $(CFLAGS) $(LIBINCLUDE) $< -o $@

clean:
	@echo -e "\033[1;31m deleting objects \033[0m" | tee make.log
	@rm -fr $(OBJECTS) $(OBJECTS_CompDensitySingle) $(OBJECTS_CompDensityDual) $(OBJECTS_IncomPressure) $(OBJECTS_CompHybrid) $(OBJECTS_PARTITION) $(OBJECTS_INITIAL) $(OBJECTS_POTENTIAL) $(OBJECTS_MapField) $(OBJECTS_MeshAMR) $(EXE) $(EXE_CompDensitySingle) $(EXE_CompDensityDual) $(EXE_IncomPressure) $(EXE_CompHybrid) $(EXE_PARTITION) $(EXE_INITIAL) $(EXE_POTENTIAL) $(EXE_MapField) $(EXE_MeshAMR) $(EXE_ExtractData) $(EXE_CompCoupled) make.log *.o
