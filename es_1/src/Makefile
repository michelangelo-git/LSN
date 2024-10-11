IDIR = ../include
CC=g++
CFLAGS=-g -I$(IDIR)

ODIR = obj
LDIR = ../lib

_DEPS = lib.h random.h random_setup.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


_OBJ = main.o lib.o random.o random_setup.o
	   
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY:run
run:
	./main

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o main
	

