/*
import sys

nodes = [1,2,3]
values = map(lambda x : 4 ** x, range(len(nodes)))

def index2assignment(ind) :
    index = ind
    tmp = [0] * len(values)
    for i in sorted(range(len(values)), reverse=True) :
        tmp[i] = index / values[i]
        index = index % values[i]
    return tmp

def assignment2index(a) :
    tmp = 0
    for i in range(len(a)) :
        tmp += (a[i] * values[i])
    return tmp


x = eval(sys.argv[1])

index = assignment2index(x)
newx = index2assignment(index)

print "%s\n%s\n%s" % (x, index, newx)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_ALLELES 4

#define TRAIT_AA 0
#define TRAIT_BA 1
#define TRAIT_AB 2
#define TRAIT_BB 3

#define NUM_NODES 4

int offsets[] = {
    1 <<  0, 
    1 <<  2, 
    1 <<  4, 
    1 <<  6, 
    1 <<  8, 
    1 << 10, 
    1 << 12,
    1 << 14,
    1 << 16,
    1 << 18
};

struct rfunction {
    int* cutset;                // eg: [1,2,3]
    int  length;
    
    float* matrix;              // len = length ** NUM_ALLELES
    
    struct rfunction* prev1;    // must be NULL if not used
    struct rfunction* prev2;    // must be NULL if not used
    
    int* cutset_indices_prev1;
    int  length_prev1;
    int* cutset_indices_prev2;
    int  length_prev2;
};

double rfunction_get(struct rfunction* rf, int index) {
    return rf->matrix[index];
}

void index2assignment(int ind, int* values, int length) {
    int index = ind;
    int i;
    
    for(i = (length - 1); i > -1; --i) {
        values[i] = index / offsets[i];
        index %= offsets[i];
    }
}

int assignment2index(int* values, int length) {
    int tmp = 0;
    int i;
    
    for(i = 0; i < length; ++i) {
        tmp += (values[i] * offsets[i]);
    }
    
    return tmp;
}

int assignment2index2(int* values, int length, int* indices, int indices_length) {
    int tmp = 0;
    int i;
    
    for(i = 0; i < indices_length; ++i) {
        tmp += (values[indices[i]] * offsets[i]);
    }
    
    return tmp;
}

int main(int argc, char** argv) {
    int values[NUM_NODES];
    int i;
    int index;
    
    if(argc != 2) {
        fprintf(stderr, "Usage: %s <index>\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    index = atoi(argv[1]);
    
    index2assignment(index, values, NUM_NODES);
    
    if(index != assignment2index(values, NUM_NODES)) {
        printf("assignment2index() does not work\n");
    }

    printf("%d = [ ", index);
    for(i = 0; i < NUM_NODES; ++i) {
        printf("%d ", values[i]);
    }
    printf("]\n");
    
    return EXIT_SUCCESS;
}

