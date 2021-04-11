all: main test


main: ./src/main.c ./src/lib.c ./src/lib.h
	gcc-10 -O3 -Wall -o ./bin/main ./src/main.c ./src/lib.c -fopenmp


test: main ./src/test.c ./src/lib.c ./src/lib.h
	gcc-10 -Wall -o ./bin/test ./src/test.c ./src/lib.c -fopenmp
	./bin/test
	python3 ./src/test.py
