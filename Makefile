.PHONY: clean
.DEFAULT_GOAL := Isaac

# executable name
EXNAME = Isaac

# compiler name
COMPILER = gfortran

# compiler flags
CFLAGS = -fopenmp                  \
         -frange-check             \
         -fbounds-check            \
         -pedantic                 \
         -Wall                     \
         -Waliasing                \
         -Wimplicit-interface      \
         -Wconversion              \
         -Wsurprising              \
         -Wunderflow               \
         -Wuninitialized           \
         -O2
#         -pg

# paths
OBJ = obj
SRC = src
BIN = bin

# objects
OBJS   = $(OBJ)/ini_reader.o       \
         $(OBJ)/global.o           \
         $(OBJ)/output.o           \
         $(OBJ)/main.o

# modules
MODULES = src/global.F90           \
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
$(OBJ)/global.o: $(SRC)/global.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/global.F90 -o $(OBJ)/global.o

$(OBJ)/ini_reader.o: $(SRC)/ini_reader.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/ini_reader.F90 -o $(OBJ)/ini_reader.o

$(OBJ)/output.o: $(SRC)/output.F90
	$(COMPILER) $(CFLAGS) -J$(OBJ) -c $(SRC)/output.F90 -o $(OBJ)/output.o

$(OBJ)/main.o: $(SRC)/main.F90
	$(COMPILER) $(CFLAGS) -I$(OBJ) -c $(SRC)/main.F90 -o $(OBJ)/main.o


# clean up
clean:
	-rm -f *~ $(EXNAME) $(OBJ)/*.o $(OBJ)/*.mod
