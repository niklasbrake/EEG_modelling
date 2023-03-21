#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char** argv){

    int N =100;
    int *neuron1_shared, *neuron2_shared;
    int *thread_shared;
    float *correlation_shared;
    int mergeCounter =0;
    int counter=0;

    int data1[N];
    int data2[N];
    for(int i = 0; i<N; i++){
        data1[i] = i;
        data2[i] = i;
    }

    int m = ceil(log2f(N+1)-1);
    int arrayLen[m];
    for(int j = 0; j<m; j++){
            arrayLen[j] = N;
            for(int i = 0; i<j; i++){
                    arrayLen[j] *= 2;
            }
    }

    // Determine how many arrays are needed to read maximal possible length of N^2
    // Initialize pointers to arrays with length 2^n
    int *neuron1_local[N], *neuron2_local[N];
    float *correlation_local[N];
    // allocate memory for first, smallest array
    neuron1_local[0] = malloc(arrayLen[0]*sizeof(int));
    neuron2_local[0] = malloc(arrayLen[0]*sizeof(int));
    correlation_local[0] = malloc(arrayLen[0]*sizeof(float));

    int arrayCounter = 0;
    int runningLength = 0;
    int totalLength = 0;
    float R;
    // Loop through all pairs of neurons
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            R = (float)i / (float)j;
            // if(R>0){
                neuron1_local[arrayCounter][runningLength] = i;
                neuron2_local[arrayCounter][runningLength] = j;
                correlation_local[arrayCounter][runningLength] = R;
                runningLength++;
                totalLength++;
                // If more memory needed, allocate next chunk of memory
                if(runningLength>=arrayLen[arrayCounter]){
                        runningLength = 0;
                        arrayCounter++;
                        neuron1_local[arrayCounter] = malloc(arrayLen[arrayCounter]*sizeof(int));
                        neuron2_local[arrayCounter] = malloc(arrayLen[arrayCounter]*sizeof(int));
                        correlation_local[arrayCounter] = malloc(arrayLen[arrayCounter]*sizeof(float));
                }
            // }
        }
    }

    // #pragma single
    // {
    // neuron1_shared = (int *) malloc(sizeof(int)*counter);
    // neuron2_shared = (int *) malloc(sizeof(int)*counter);
    // thread_shared = (int *) malloc(sizeof(int)*counter);
    // correlation_shared = (float *) malloc(sizeof(int)*counter);
    // }

    // #pragma omp critical
    // {
    //     int M;
    //     for(int j = 0; j<=arrayCounter; j++){
    //         if(j<arrayCounter){
    //             M = arrayLen[j];
    //         }
    //         else{
    //             M = runningLength;
    //         }
    //         for(int i = 0; i<M; i++){
    //             neuron1_shared[mergeCounter] = neuron1_local[j][i];
    //             neuron2_shared[mergeCounter] = neuron2_local[j][i];
    //             correlation_shared[mergeCounter] = correlation_local[j][i];
    //             thread_shared[mergeCounter] = omp_get_thread_num();
    //             mergeCounter++;
    //         }
    //     }
    // }
    for(int j = 0; j<m; j++){
        printf("%d\n",arrayLen[j]);
    }
    printf("%d\n",counter);
    // free(neuron1_local);
    // free(neuron2_local);
    // free(correlation_local);
    // for(int i = 0; i < mergeCounter; i++){
    //     printf("%d, %d, %.4f (thread %d)\n", neuron1_shared[i],neuron2_shared[i],correlation_shared[i],thread_shared[i]);
    // }
    return 0;
}