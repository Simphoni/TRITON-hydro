CC := g++

EXECS= bin2ascii_rmap

all: $(EXECS)

bin2ascii:
	$(CC) -o $@ bin2ascii_rmap.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ bin2ascii_rmap
