f90files = $(wildcard *.f90)
OBJ= $(f90files:.f90=.o)

F90 = gfortran  -g 
F77 = gfortran -g 
LINK = $(F77)
# CFLAGS = -cpp -fsloppy-char
PFLAGS =  -cpp $(CFLAGS)

LIBS = /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf

EXEC = arttest.exe
all:$(EXEC)

$(EXEC): $(OBJ)
	$(LINK) $(CFlAGS) $(OBJ) -o $(EXEC)  $(LIBS)


%.o: %.f90
	        $(F90) $(PFLAGS) -c  -o $@ $*.f90

# Dependencies
art.o : art.f90 defs.o random.o art_lanczos.o
art_step.o: art_step.f90  find_saddle.o
find_saddle.o : find_saddle.f90

clean:
	rm *.o *.mod $(EXEC)
