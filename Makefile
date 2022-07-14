### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE         #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0    #Positions in Kpc/h
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DSTORE_IDS
EXTRAS += -DCOMPUTE_EP
#EXTRAS += -DTYPE_PART=0
#EXTRAS += -DTYPE_TWO_GADGET
#EXTRAS += -DFILE_ASCII
#EXTRAS += -DCHANGE_POSITION
#EXTRAS += -DWRITE_WEIGHT
EXTRAS += -DPRUNED
EXTRAS += -DBRANCH_SURVIVE
EXTRAS += -DSORT
#EXTRAS += -DLEVEL_PRUNED=10
EXTRAS += -DITERA

#CXX
CXX     := g++ 
DC     := -DNTHREADS=16
DC     += -DLOCK
GSLL   := -lgsl -lgslcblas
CFLAGS := -Wall -O3 -march=native -ftree-vectorize -ansi -pedantic -g -fopenmp
LIBS   := -lm $(GSLL) 
VPP_INC := -I./voro++
VPP_LIB := -L./voro++ -lvoro++
#VPP_INC := -I/home/lpereyra/voro++/include/voro++
#VPP_LIB := -L/home/lpereyra/voro++/lib -lvoro++
LDFLAGS = -lm -lgfortran 
FC      = gfortran
FFLAGS  = -cpp -O3 -fbounds-check 

.PHONY : clean cleanall

MAKEFILE := Makefile

OBJS :=  variables.o io.o allocate.o leesnap.o grid.o \
iden.o properties.o deltas.o octree.o \
voronoi.o kruskal.o select.o libfitpack.o vol.o properties_fil.o 

OBJS_FITPACK := ./fitpack/parcur.o \
							  ./fitpack/fpback.o \
						    ./fitpack/fpbspl.o \
						    ./fitpack/fpchec.o \
						    ./fitpack/fppara.o \
						    ./fitpack/fpdisc.o \
						    ./fitpack/fpgivs.o \
						    ./fitpack/fpknot.o \
						    ./fitpack/fprati.o \
						    ./fitpack/fprota.o \
							  ./fitpack/splev.o  

HEADERS := $(patsubst %.o,$(OBJS),$(OBJS_FITPACK))

EXEC := main.x

todo: $(EXEC)

%.f.o: $(OBJS_FITPACK) 
	$(FC) $(FFLAGS) -c $<

%.o: %.cc %.hh $(MAKEFILE)
	$(CXX) $(EXTRAS) $(VPP_INC) $(CFLAGS) $(DC) -c $<

main.x: main.cc $(OBJS) $(OBJS_FITPACK)
	$(MAKE) -C voro++
	$(CXX) $(CFLAGS) $(EXTRAS) $^ -o $@ $(LIBS) $(VPP_LIB) $(LDFLAGS)

clean:
	rm -rf $(OBJS_FITPACK)
	rm -rf $(OBJS)
	rm -rf main.o

cleanall:
	$(MAKE) -C voro++ clean
	rm -rf $(OBJS_FITPACK)
	rm -rf $(OBJS)
	rm -rf main.o
	rm -rf $(EXEC)
