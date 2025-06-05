.PHONY: clean
.DEFAULT_GOAL := Isaac

# executable name
EXNAME = Isaac

# compiler name
#COMPILER = /usr/lib64/mvapich2/bin/mpif90
#COMPILER = mpif90
COMPILER = gfortran

# http://gcc.gnu.org/onlinedocs/gfortran/Fortran-Dialect-Options.html
# compiler flags
CFLAGS = -fopenmp                  \
         -ffree-line-length-none   \
         -frange-check             \
         -fbounds-check            \
         -funroll-loops            \
         -Wall                     \
         -Wextra                   \
         -Waliasing                \
         -Wconversion              \
         -Wsurprising              \
         -Wunderflow               \
         -Wuninitialized           \
         -cpp                      \
         -O3
#         -Wimplicit-interface                      
#         -pedantic
#         -g
#         -pg

# paths
OBJ = obj
SRC = src
BIN = bin

# objects
OBJS   = $(OBJ)/precision.o        \
         $(OBJ)/strings.o          \
         $(OBJ)/ini_file_reader.o  \
         $(OBJ)/global.o           \
         $(OBJ)/output.o           \
         $(OBJ)/input.o            \
         $(OBJ)/initialize.o       \
         $(OBJ)/integrators.o      \
         $(OBJ)/main.o

# modules
MODULES = src/precision.F90        \
          src/strings.F90          \
          src/global.F90           \
          src/output.F90           \
          src/input.F90            \
          src/initialize.F90       \
          src/integrators.F90      \
          src/ini_reader.F90

# create the executable 
Isaac: make_dirs $(OBJS)
	$(COMPILER) $(CFLAGS) -o $(EXNAME) $(OBJS)

# create a new executable 
new: clean Isaac

# make dirs
make_dirs:
	mkdir -p obj

# compile the source files
$(OBJ)/precision.o: $(SRC)/precision.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/precision.F90 -o $(OBJ)/precision.o

$(OBJ)/strings.o: $(SRC)/strings.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/strings.F90 -o $(OBJ)/strings.o

$(OBJ)/ini_file_reader.o: $(SRC)/ini_file_reader.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/ini_file_reader.F90 -o $(OBJ)/ini_file_reader.o

$(OBJ)/global.o: $(SRC)/global.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/global.F90 -o $(OBJ)/global.o

$(OBJ)/output.o: $(SRC)/output.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/output.F90 -o $(OBJ)/output.o

$(OBJ)/input.o: $(SRC)/input.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/input.F90 -o $(OBJ)/input.o

$(OBJ)/initialize.o: $(SRC)/initialize.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/initialize.F90 -o $(OBJ)/initialize.o

$(OBJ)/integrators.o: $(SRC)/integrators.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/integrators.F90 -o $(OBJ)/integrators.o

$(OBJ)/main.o: $(SRC)/main.F90
	$(COMPILER) $(CFLAGS) -I$(OBJ) -c $(SRC)/main.F90 -o $(OBJ)/main.o


# clean up
clean:
	-rm *~ $(EXNAME) $(OBJ)/*.o $(OBJ)/*.mod
