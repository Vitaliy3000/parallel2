#include "parser.h"


char* MAX_SIZE_OF_STRING_ERROR = "Max size of argumet %s must be less 40 letters\n";
char* DUPLICATE_ARGUMENT_ERROR = "Duplicate argument of %s\n";
char* NOT_FOUND_ARGUMENT_ERROR = "Not found argument of %s\n";

char* LOW_ARGUMENT_ERROR_INTEGER = "Argument %s must be bigger %d\n";
char* HIGH_ARGUMENT_ERROR_INTEGER = "Argument %s must be lower %d\n";
char* LOW_ARGUMENT_ERROR_FLOAT = "Argument %s must be bigger %.2f\n";
char* HIGH_ARGUMENT_ERROR_FLOAT = "Argument %s must be lower %.2f\n";
char* INTEGER_ARGUMENT_ERROR = "Argument %s must be ineger\n";
char* FLOAT_ARGUMENT_ERROR = "Argument %s must be double\n";


int strict_atoi(char str[], int *ptr) {
    for(int i = 0; str[i] != '\0'; i++) {
        if (((str[i] > '9') || (str[i] < '0')) && (str[i] != '-')) {
            return -1;
        }
    }

    *ptr = atoi(str);
    return 0;
}


double strict_atof(char str[], double *ptr) {
    for(int i = 0; str[i] != '\0'; i++) {
        if (((str[i] > '9') || (str[i] < '0')) && (str[i] != '-') && (str[i] != '.') && (str[i] != 'e')) {
            return -1;
        }
    }

    *ptr = atof(str);
    return 0;
}


void parse_int_arg(int argc, char *argv[], char expected_key[], int ptr[], int low, int high) {
    char value[MAX_SIZE_OF_STRING] = "";
    
    for (int i = 1; i < argc-1; i++) {
        if (strcmp(argv[i], expected_key) == 0) {
            if (value[0] != '\0') {
                fprintf(stderr, DUPLICATE_ARGUMENT_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            if (strlen(argv[i+1]) >= MAX_SIZE_OF_STRING) {
                fprintf(stderr, MAX_SIZE_OF_STRING_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            strcpy(value, argv[i+1]);
        }
    }

    if (value[0] == '\0') {
        fprintf(stderr, NOT_FOUND_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (strict_atoi(value, ptr) == -1) {
        fprintf(stderr, INTEGER_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr < low) {
        fprintf(stderr, LOW_ARGUMENT_ERROR_INTEGER, expected_key, low);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr > high) {
        fprintf(stderr, HIGH_ARGUMENT_ERROR_INTEGER, expected_key, high);
        ERROR_FLAG = 1;
        return;
    }

    return;
}


void parse_double_arg(int argc, char *argv[], char expected_key[], double ptr[], double low, double high) {
    char value[MAX_SIZE_OF_STRING] = "";

    for (int i = 1; i < argc-1; i++) {
        if (strcmp(argv[i], expected_key) == 0) {
            if (value[0] != '\0') {
                fprintf(stderr, DUPLICATE_ARGUMENT_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            if (strlen(argv[i+1]) >= MAX_SIZE_OF_STRING) {
                fprintf(stderr, MAX_SIZE_OF_STRING_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            strcpy(value, argv[i+1]);
        }
    }

    if (value[0] == '\0') {
        fprintf(stderr, NOT_FOUND_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (strict_atof(value, ptr) == -1) {
        fprintf(stderr, FLOAT_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr < low) {
        fprintf(stderr, LOW_ARGUMENT_ERROR_FLOAT, expected_key, low);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr > high) {
        fprintf(stderr, HIGH_ARGUMENT_ERROR_FLOAT, expected_key, high);
        ERROR_FLAG = 1;
        return;
    }

    return;
}


void parse_str_arg(int argc, char *argv[], char expected_key[], char *str) {
    for (int i = 1; i < argc-1; i++) {
        if (strcmp(argv[i], expected_key) == 0) {
            if (str[0] != '\0') {
                fprintf(stderr, DUPLICATE_ARGUMENT_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            if (strlen(argv[i+1]) >= MAX_SIZE_OF_STRING) {
                fprintf(stderr, MAX_SIZE_OF_STRING_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            strcpy(str, argv[i+1]);
        }
    }
}