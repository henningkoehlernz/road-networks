CC = g++ -std=c++20 -O3 -Wall -Wextra -o
TCC = g++ -std=c++20 -ggdb -Wall -Wextra -o
INC = src/road_network.cpp src/util.cpp

default: main
all: main test
main:
	$(CC) cut src/main.cpp $(INC)
test:
	$(TCC) test src/test.cpp $(INC)
clean:
	rm cut test

.PHONY: main test
