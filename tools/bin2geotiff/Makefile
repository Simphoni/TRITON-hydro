CC := g++

EXECS= bin2geotiff

all: $(EXECS)

bin2geotiff:
	$(CC) -o $@ bin2geotiff.cpp -O3 --std=c++11 -lgdal -Wunused-result

.PHONY: clean
clean:
	rm -f *.so *.o *~ bin2geotiff
