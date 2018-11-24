#include <stdio.h>

#define ACONST 10
#define VECSIZE 4

int fmult(int a, int b);
int r1;
int vector[VECSIZE]={1,2,3,4}; // This declares a vector of size VECSIZE.

int main(){
    int a=1, b=2, c; // Declaration of local variables
    c = a + b + ACONST;
    r1 = fmult(12, 13);
    printf("a + b + ACONST is %d \nr1 = %d \n", c , r1);
    getchar();
    return(0); // Return 0
}

int fmult(int a, int b){
    int c;
    int i;
    c = a * b;
    for (i=0; i<VECSIZE; i++){
        c += vector[i]*vector[i];
    };
    return(c);
}