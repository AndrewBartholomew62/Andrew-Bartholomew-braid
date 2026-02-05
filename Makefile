
OBJFILES = util.o debug.o bigint.o class-control.o input.o preprocessor.o\
           braid-util.o gauss-to-peer.o generic-code-io.o generic-code-util.o gauss-orientation.o reidemeister.o \
           main.o bracket.o braid.o  braidfns.o generic-code.o hamiltonian.o vogel.o vogelfns.o homology.o                    

PREP_OBJFILES = debug.o preprocessor-main.o preprocessor.o          

TEST_OBJFILES = debug.o run-test.o preprocessor.o          
         
DEPS     = ./include/* 

COMPILE  = g++ -Wall -Wno-misleading-indentation -std=c++11 -I ./include -g -c $< -o $@ 

all: $(OBJFILES)
	g++ -o braid $(OBJFILES) 

preprocess: $(PREP_OBJFILES)
	g++ -o preprocess $(PREP_OBJFILES)

run-test: $(TEST_OBJFILES)
	g++ -o run-test $(TEST_OBJFILES)

util.o: ./src/util.cpp $(DEPS)
	$(COMPILE)

debug.o: ./src/debug.cpp $(DEPS)
	$(COMPILE)

bigint.o: ./src/bigint.cpp ./include/bigint-scalar.h
	$(COMPILE)

class-control.o: ./src/class-control.cpp $(DEPS)
	$(COMPILE)

input.o: ./src/input.cpp $(DEPS)
	$(COMPILE)

braid-util.o: ./src/braid-util.cpp $(DEPS)
	$(COMPILE)

gauss-to-peer.o: ./src/gauss-to-peer.cpp $(DEPS)
	$(COMPILE)

generic-code-io.o: ./src/generic-code-io.cpp $(DEPS)
	$(COMPILE)

generic-code-util.o: ./src/generic-code-util.cpp $(DEPS)
	$(COMPILE)

gauss-orientation.o: ./src/gauss-orientation.cpp $(DEPS)
	$(COMPILE)
	
hamiltonian.o: ./src/hamiltonian.cpp $(DEPS)
	$(COMPILE)	

reidemeister.o: ./src/reidemeister.cpp $(DEPS)
	$(COMPILE)

main.o: ./src/main.cpp $(DEPS)
	$(COMPILE)

braid.o: ./src/braid.cpp $(DEPS)
	$(COMPILE)

braidfns.o: ./src/braidfns.cpp $(DEPS)
	$(COMPILE)

generic-code.o: ./src/generic-code.cpp $(DEPS)
	$(COMPILE)

bracket.o: ./src/bracket.cpp $(DEPS)
	$(COMPILE)


vogel.o: ./src/vogel.cpp $(DEPS)
	$(COMPILE)

vogelfns.o: ./src/vogelfns.cpp $(DEPS)
	$(COMPILE)

homology.o: ./src/homology.cpp $(DEPS)
	$(COMPILE)

preprocessor.o: ./src/preprocessor.cpp $(DEPS)
	$(COMPILE)

preprocessor-main.o: ./src/preprocessor-main.cpp $(DEPS)
	$(COMPILE)

run-test.o: ./src/run-test.cpp $(DEPS)
	$(COMPILE)
	
clean:
	rm -f $(OBJFILES)
