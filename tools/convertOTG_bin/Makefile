CC := g++

EXECS= headerBinary

all: $(EXECS)

bin2ascii:
	$(CC) -o $@ headerBinary.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ headerBinary
