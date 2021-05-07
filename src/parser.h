#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// #ifndef PARSER_H
// #define PARSER_H

#define MAX_SIZE_OF_STRING 40


static int ERROR_FLAG = 0;

int strict_atoi(char[], int*);
double strict_atof(char[], double*);
void parse_int_arg(int, char*[], char[], int[], int, int);
void parse_double_arg(int, char*[], char[], double[], double, double);
void parse_str_arg(int, char*[], char[], char*);

// #endif