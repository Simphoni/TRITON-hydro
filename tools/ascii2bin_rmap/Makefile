CC := g++

EXECS= ascii2bin_rmap

all: $(EXECS)

ascii2bin:
	$(CC) -o $@ ascii2bin_rmap.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ ascii2bin_rmap
