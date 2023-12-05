CC := g++

EXECS= bin2ascii_dem

all: $(EXECS)

bin2ascii_dem:
	$(CC) -o $@ bin2ascii_dem.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ bin2ascii_dem