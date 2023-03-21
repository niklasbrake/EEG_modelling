#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char** argv){
        N = 100;
        int *pointers[5];
        pointers[0] = malloc(N*sizeof(int));
        for(i=0; i<N; i++){
                pointers[0][i] = i;
        }
        for(j=0; j<5; j++){
                for(i=0; i<N; i++){
                        printf("%d\n",pointers[j][i]);
                }
        }
        return 0;
}