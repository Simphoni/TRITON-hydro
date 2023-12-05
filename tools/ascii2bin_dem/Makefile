CC := g++

EXECS= ascii2bin_dem

all: $(EXECS)

ascii2bin_dem:
	$(CC) -o $@ ascii2bin_dem.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ ascii2bin_dem