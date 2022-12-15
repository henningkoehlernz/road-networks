CC = g++ -O3 -Wall -o
INC = src/road_network.cpp

default: main
all: main
main:
	$(CC) cut src/cut_index.cpp $(INC)
clean:
	rm cut

.PHONY: main
