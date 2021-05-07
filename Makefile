all: task1 task2 test


task1: ./src/task1.c ./src/matrix.h ./src/matrix.c ./src/lib.h ./src/lib.c ./src/parser.h ./src/parser.c ./src/kernals_omp.h ./src/kernals_omp.c
	gcc-10 -O3 -Wall -fopenmp -o ./bin/task1 ./src/task1.c ./src/matrix.c ./src/lib.c ./src/parser.c ./src/kernals_omp.c


task2: ./src/task2.c ./src/matrix.h ./src/matrix.c ./src/lib.h ./src/lib.c ./src/parser.h ./src/parser.c ./src/kernals.h ./src/kernals.c
	mpicc -O3 -Wall -o ./bin/task2 ./src/task2.c ./src/matrix.c ./src/lib.c ./src/parser.c ./src/kernals.c


test: task1 task2 ./src/test_kernals.c ./src/test_kernals_omp.c
	gcc-10 -Wall -o ./bin/test_kernals ./src/test_kernals.c ./src/kernals.c -fopenmp
	gcc-10 -Wall -o ./bin/test_kernals_omp ./src/test_kernals_omp.c ./src/kernals_omp.c -fopenmp
	./bin/test_kernals
	./bin/test_kernals_omp
	python3 ./src/test_tasks.py
