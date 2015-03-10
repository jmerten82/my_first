#Choose code package you want to compile

#SELECTION	= "MPIField"
#SELECTION	= "singleCPUHighres"
#SELECTION	= "MPIReconstruction"
#SELECTION	= "MPIFieldBS"
#SELECTION 	= "MPIReconstructionBS"
SELECTION	= "singleCPUHighresBS"




ifeq ($(SELECTION),"MPIField")
NAME	 = SAW_FIELD
CC 	 = mpiCC
CCFLAGS	 =  
LDFLAGS	 = -lclapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -lCCfits
ADDLIB	 = -L/usr/lib64/atlas/
ADDINC	 = -I./ -I/usr/lib64/atlas/ -I/usr/include/  		
DRIVER   = sawmpi_griddriver.cpp
SRC      = mpi_inputfield.cpp rw_fits.cpp options.cpp soph_math.cpp littlelibastro.cpp util.cpp masked_fin_dif.cpp mpi_comm.cpp
endif


ifeq ($(SELECTION),"MPIReconstruction")
NAME	 = SAW_REC
CC	 = mpiCC
CCFLAGS	 = 
LDFLAGS	 = -lclapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -lCCfits
ADDLIB   = -L/usr/lib64/atlas/
ADDINC   = -I./ -I/usr/lib64/atlas/ -I/usr/include/
DRIVER	 = sawmpi_recdriver.cpp
SRC	 = galaxycluster.cpp options.cpp masked_fin_dif.cpp mpi_comm.cpp mpi_coreroutines.cpp mpi_datacontainer.cpp rw_fits.cpp interpol.cpp littlelibastro.cpp util.cpp soph_math.cpp
endif

ifeq ($(SELECTION),"singleCPUHighres")
NAME	 = SAW_HD
CC	 = g++
CCFLAGS	 =
LDFLAGS	 = -lclapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -lCCfits
ADDLIB   = -L/usr/lib64/atlas/
ADDINC   = -I./ -I/usr/lib64/atlas/ -I/usr/include/
DRIVER   = saw_highres.cpp
SRC	 = masked_fin_dif.cpp prep_strong.cpp rw_fits.cpp littlelibastro.cpp util.cpp reconstruction_kernel.cpp galaxycluster.cpp options.cpp interpol.cpp soph_math.cpp 
endif

ifeq ($(SELECTION),"MPIFieldBS")
NAME	 = SAW_FIELD_BS
CC 	 = mpiCC
CCFLAGS	 =  
LDFLAGS	 = -lclapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -lCCfits 
ADDLIB	 = -L/usr/lib64/atlas/
ADDINC	 = -I./ -I/usr/lib64/atlas/ -I/usr/include/  		
DRIVER   = sawmpi_griddriverBS.cpp
SRC      = mpi_inputfield.cpp rw_fits.cpp options.cpp soph_math.cpp littlelibastro.cpp util.cpp masked_fin_dif.cpp mpi_comm.cpp bootstrap.cpp cat_reader.cpp cat_writer.cpp
endif

ifeq ($(SELECTION),"MPIReconstructionBS")
NAME	 = SAW_REC_BS
CC	 = mpiCC
CCFLAGS	 = 
LDFLAGS	 = -lclapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -lCCfits 
ADDLIB   = -L/usr/lib64/atlas/
ADDINC   = -I./ -I/usr/lib64/atlas/ -I/usr/include/ 
DRIVER	 = sawmpi_recdriverBS.cpp
SRC	 = galaxycluster.cpp options.cpp masked_fin_dif.cpp mpi_comm.cpp mpi_coreroutines.cpp mpi_datacontainer.cpp rw_fits.cpp interpol.cpp littlelibastro.cpp util.cpp soph_math.cpp 
endif

ifeq ($(SELECTION),"singleCPUHighresBS")
NAME	 = SAW_HD_BS
CC	 = g++
CCFLAGS	 =
LDFLAGS	 = -lclapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -lCCfits
ADDLIB   = -L/usr/lib64/atlas
ADDINC   = -I./ -I/usr/lib64/atlas/ -I/usr/include/  
DRIVER   = saw_highresBS.cpp
SRC	 = masked_fin_dif.cpp prep_strong.cpp rw_fits.cpp littlelibastro.cpp util.cpp reconstruction_kernel.cpp galaxycluster.cpp options.cpp interpol.cpp soph_math.cpp
endif



HEADER	= $(SRC:.cpp=.h)
OBJECT	= $(SRC:.cpp=.o)
NR	= nr.h nrexit.cpp nrtypes.h nrtypes_lib.h nrtypes_nr.h nrutil.h nrutil_tnt.h nrutil_val.h 
INCLUDES = $(GSLINC) $(FITSINC) $(ATLASINC) $(ADDINC)
LIBS = $(GSLLIB) $(FITSLIB) $(ATLASLIB) $(ADDLIB)

.PHONY: clean uninstall clearance

%.o: %.cpp %.h
	$(CC) $(CCFLAGS) -c $< $(INCLUDES) $(LIBS)
%.o: %.cpp
	$(CC) $(CCFLAGS) -c $< $(INCLUDES) $(LIBS)

$(NAME): $(DRIVER:.cpp=.o) $(OBJECT) $(HEADER)
	$(CC) $(CCFFLAGS) -o $@ $< $(OBJECT) $(INCLUDES) $(LIBS) $(LDFLAGS)

default: $(NAME)

all:	 $(NAME) clearance

clean:
	   mkdir src; \
           mv $(SRC) $(HEADER) $(DRIVER) ./src; \
	   mv $(NR) ./src;\
	   rm *.cpp; \
	   rm *.h;\
	   rm -rf $(OBJECT) $(NAME) $(DRIVER:.cpp=.o)

install: $(NAME)

uninstall: clean
	   rm -rf ./src; \


saw_griddriver.o: options.o inputfield.o
sawmpi_recdriver.o: options.o rw_fits.o mpi_comm.o mpi_coreroutines.o mpi_datacontainer.o interpol.o soph_math.o util.o
galaxycluster.o: util.o rw_fits.o options.o
inputfield.o: rw_fits.o util.o soph_math.o options.o masked_fin_dif.o
prep_strong.o: interpol.o options.o util.o
reconstruction_kernel.o: galaxycluster.o masked_fin_dif.o options.o util.o
util.o: littlelibastro.o rw_fits.o
masked_fin_dif.o: util.o
interpol.o: util.o
mpi_inputfield.o: rw_fits.o util.o soph_math.o options.o masked_fin_dif.o mpi_comm.o
mpi_coreroutines.o: options.o galaxycluster.o util.o rw_fits.o mpi_comm.o mpi_datacontainer.o
mpi_datacontainer.o: options.o rw_fits.o mpi_comm.o masked_fin_dif.o util.o
bootstrap.o: cat_writer.o
cat_writer.o: cat_reader.o
cat_reader.o: util.o
