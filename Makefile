all: main test


main: ./src/main.c ./src/lib.c ./src/lib.h
	gcc -O -Wall -o ./bin/main ./src/main.c ./src/lib.c -fopenmp


test: main ./src/test.c ./src/lib.c ./src/lib.h
	gcc -Wall -o ./bin/test ./src/test.c ./src/lib.c -fopenmp
	./bin/test.exe
	python ./src/test.py
