
OBJFILES = util.o debug.o bigint-scalar.o class-control.o \
           gauss-orientation.o reidemeister.o \
           main.o braid.o generic-code.o bracket.o braidfns.o vogel.o vogelfns.o
         
DEPS     = ./include/* 

COMPILE  = g++ -Wall -Wno-misleading-indentation -std=c++11 -I ./include -g -c $< -o $@ 

all: $(OBJFILES)
	g++ -o braid $(OBJFILES) 

util.o: ./src/util.cpp $(DEPS)
	$(COMPILE)

debug.o: ./src/debug.cpp $(DEPS)
	$(COMPILE)

bigint-scalar.o: ./src/bigint-scalar.cpp $(DEPS)
	$(COMPILE)

class-control.o: ./src/class-control.cpp $(DEPS)
	$(COMPILE)

gauss-orientation.o: ./src/gauss-orientation.cpp $(DEPS)
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

clean:
	rm -f $(OBJFILES)
