CC := g++

EXECS= discharge2velocity

all: $(EXECS)

discharge2velocity:
	$(CC) -o $@ discharge2velocity.cpp -g --std=c++11

.PHONY: clean
clean:
	rm -f *.so *.o *~ discharge2velocity