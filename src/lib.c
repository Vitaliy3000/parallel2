#include "lib.h"


void* strict_malloc(size_t size) {
    void *new_ptr = malloc( size);
 
    if (new_ptr == NULL) {
        perror("malloc return NULL");
        exit(EXIT_FAILURE);
    }

    return new_ptr;
}


void* strict_calloc(size_t num, size_t size) {
    void *new_ptr = calloc(num, size);
 
    if (new_ptr == NULL) {
        perror("calloc return NULL");
        exit(EXIT_FAILURE);
    }

    return new_ptr;
}


void* strict_realloc(void *ptr, size_t newsize) {
    void *new_ptr = realloc(ptr, newsize);

    if (new_ptr == NULL) {
        perror("realloc return NULL");
        exit(EXIT_FAILURE);
    }

    return new_ptr;
}
