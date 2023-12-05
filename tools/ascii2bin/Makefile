CC := g++

EXECS= ascii2bin

all: $(EXECS)

ascii2bin:
	$(CC) -o $@ ascii2bin.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ ascii2bin