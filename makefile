CC = g++ -std=c++20 -O3 -Wall -Wextra -o
TCC = g++ -std=c++20 -ggdb -Wall -Wextra -o
INC = src/road_network.cpp src/util.cpp

default: main
all: main index query topcut test
main:
	$(CC) cut src/main.cpp $(INC)
index:
	$(CC) index src/index.cpp $(INC)
query:
	$(CC) query src/query.cpp $(INC)
topcut:
	$(CC) topcut src/topcut.cpp $(INC)
test:
	$(TCC) test src/test.cpp $(INC)
clean:
	rm cut index query topcut test

.PHONY: main index query topcut test
