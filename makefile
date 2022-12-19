CC = g++ -std=c++20 -O3 -Wall -o
INC = src/road_network.cpp

default: main
all: main test
main:
	$(CC) cut src/cut_index.cpp $(INC)
test:
	$(CC) test src/test.cpp $(INC)
clean:
	rm cut test

.PHONY: main test
