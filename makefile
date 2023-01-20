CC = g++ -std=c++20 -O3 -Wall -Wextra -o
TCC = g++ -std=c++20 -ggdb -Wall -Wextra -o
INC = src/road_network.cpp src/util.cpp

default: main
all: main topcut test
main:
	$(CC) cut src/main.cpp $(INC)
topcut:
	$(CC) topcut src/topcut.cpp $(INC)
test:
	$(TCC) test src/test.cpp $(INC)
clean:
	rm cut test topcut

.PHONY: main topcut test
