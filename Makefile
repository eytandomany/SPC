# Makefile to clustering using SW 

CC = gcc

CFLAGS =  -O3 
EFILE =  ./SW

OBJ  =  SW.o RaggedArray.o aux1.o aux2.o edge.o distance.o io.o param.o timer.o utilities.o

all: rmexec $(EFILE)

tmp: $(OBJ) tmp.o
	$(CC) utilities.o tmp.o $(CFLAGS) -o tmp -lm

$(EFILE): $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o $(EFILE) -lm

SW.o aux1.o aux2.o edge.o distance.o io.o: SW.h
timer.o: timer.h
param.o: param.h
RaggedArray.o: RaggedArray.h
utilities.o: utilities.h
SW.h: timer.h param.h RaggedArray.h utilities.h
	touch SW.h

rmexec:
	rm -f $(EFILE)

clean: 
	@echo "Removing object files ..."
	@rm $(OBJ)
	@echo "Done"







