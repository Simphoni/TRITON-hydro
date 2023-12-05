CC := g++

EXECS= bin2ascii

all: $(EXECS)

bin2ascii:
	$(CC) -o $@ bin2ascii.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ bin2ascii