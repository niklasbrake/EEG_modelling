#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int computeFileSize(char *filename, int* saveWidth, int* saveHeight)
{
    FILE* stream = fopen(filename, "r");
    if(!stream){
        printf("ERROR: %s\n",filename);
    }
    char line[100000];
    int nWidth = 0;
    int nHeight = 0;
    char *token;
    while (fgets(line, 100000, stream))
    {
        token = strtok(line, ",");
        int sum = 0;
        while( token != NULL )
        {
            sum++;
            token = strtok(NULL, ",");
        }
        if(sum>nWidth)
        {
            nWidth = sum;
        }
        nHeight++;
    }
    fclose(stream);
    *saveWidth = nWidth;
    *saveHeight = nHeight-1;
    return 0;
}

int *readSpikeTrains(char *filename,int nWidth,int nHeight)
{
    int *data;
    data = malloc(sizeof(int) * nWidth*nHeight*5);
    if (!data) {
        perror("Error allocating memory");
        abort();
    }
    memset(data, 0, sizeof(int)*nWidth*nHeight*5);

    FILE* stream = fopen(filename, "r");
    FILE* errorFile = fopen("./errors.txt","w");
    char line[100000];
    int i = 0;
    int j = 0;
    while (fgets(line, 100000, stream))
    // for (int i = 0; i<100; i++)
    {
        char *token;
        // First token should be e/i so ignore
        token = strtok(line, ",");
        token = strtok(NULL, ",");
        j=0;
        while( token != NULL )
        {
            // First token should be e/i so ignore
            data[nWidth*i+j] = atoi(token);
            token = strtok(NULL, ",");
            j++;
        }
        i++;
    }
    // fprintf(errorFile,"%d,%d\n",sizeof(int) * nWidth*nHeight,nWidth*i+j);
    fclose(stream);
    fclose(errorFile);

    return data;
}

int computePairwiseCorrelation(int *data,int nWidth, int nHeight, char*saveFileName)
{
    int dt = 10;
    int m = ceil(log2f(nHeight+1)-1);
    int arrayLen[m];
    for(int j = 0; j<m; j++){
        arrayLen[j] = nHeight;
        for(int i = 0; i<j; i++){
            arrayLen[j] *= 2;
        }
    }


    FILE* outputFile = fopen(saveFileName,"w");

    // Estimate length of spike trains
    float totalTime = 0;
    int i0;
    for (int i = 0; i<nHeight; i++) {
        i0 = i*nWidth;
        while(data[i0]!=0) {
            if(data[i0]>totalTime){
                totalTime = data[i0];
            }
            i0++;
        }
    }

    // TO DO: implement with balanced binary tree to make more efficient...
    float T[nHeight];
    for (int i = 0; i<nHeight; i++) {
        // Compute T
        i0 = i*nWidth;

        T[i] = 0;
        int t0 = 0;
        int t1 = 0;
        while(data[i0]!=0) {
            t0 = data[i0]-dt;
            if(t0<t1){
                // Shorten current interval to abut previous
                t0 = t1;
            }
            t1 = data[i0]+dt;
            if(t1>totalTime){
                t1 = totalTime;
            }
            T[i] = T[i]+t1-t0;
            i0++;
        }
        T[i] = T[i]/totalTime;
    }


    #pragma omp parallel
    {
        // Initialize pointers to arrays with length 2^n
        int *neuron1_local[m], *neuron2_local[m];
        int *data_local = data;
        float *T_local = T;
        float *correlation_local[m];
        // allocate memory for first, smallest array
        neuron1_local[0] = malloc(arrayLen[0]*sizeof(int));
        neuron2_local[0] = malloc(arrayLen[0]*sizeof(int));
        correlation_local[0] = malloc(arrayLen[0]*sizeof(float));

        float R;
        int i;
        int j;
        int i0;
        int j0;
        int overlapCounter;
        int overlapCounter2;
        int totalCounter;
        int arrayCounter = 0;
        int runningLength = 0;
        int totalLength = 0;
        int *j0Counter[nWidth];
        int k;
        int totalJspikes = 0;
        float STTC;

        #pragma omp for
        for (i = 0; i<nHeight; i++) {
            for (j = i+1; j<nHeight; j++) {
                i0 = i*nWidth;
                j0 = j*nWidth;

                totalCounter = 0;
                int Pi = 0;
                int Pj = 0;

                memset(j0Counter, 0, sizeof(int)*nWidth);
                while(data_local[i0]!=0) {
                    j0 = j*nWidth;
                    overlapCounter = 0;
                    k = 0;
                    while(data_local[j0]!=0) {
                        totalCounter++;
                        if(abs(data_local[i0]-data_local[j0])<=8) {
                            // How many overlaps with spike j0?
                            j0Counter[k]++;
                            overlapCounter++;
                        }
                        j0++;
                        k++;
                    }
                    // If spike i0 overlapped with at least 1 j0
                    if(overlapCounter>0){
                        Pi++;
                    }
                    i0++;
                }
                totalJspikes = (j0-j*nWidth);
                for(int k1=0; k1<totalJspikes; k1++){
                    if(j0Counter[k1]>0){
                        Pj++;
                    }
                }
                if(totalCounter>0){
                    float Pi2 = (float)Pi / (float)(i0-i*nWidth);
                    float Pj2 = (float)Pj / (float)totalJspikes;
                    STTC = (Pi2-T_local[j])/(1-Pi2*T_local[j]) + (Pj2-T_local[i])/(1-Pj2*T_local[i]);
                    if(STTC>0.1){
                    // if(STTC>0){
                        neuron1_local[arrayCounter][runningLength] = i;
                        neuron2_local[arrayCounter][runningLength] = j;
                        correlation_local[arrayCounter][runningLength] = STTC/2;
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
                    }
                }
            }
        }
        #pragma omp critical
        {
            int M;
            for(int j = 0; j<=arrayCounter; j++){
                if(j<arrayCounter){
                    M = arrayLen[j];
                }
                else{
                    M = runningLength;
                }
                for(int i = 0; i<M; i++){
                    fprintf(outputFile,"%d,%d,%.6f\n",neuron1_local[j][i],neuron2_local[j][i],correlation_local[j][i]);
                }
            }
        }
    }
    fclose(outputFile);
    return 0;
}

int main(int argc, char **argv)
{
    int nWidth,nHeight;

    char *inputFile = argv[1];
    char *saveFileName = argv[2];

    computeFileSize(inputFile,&nWidth,&nHeight);
    int *data = readSpikeTrains(inputFile,nWidth,nHeight);
    computePairwiseCorrelation(data,nWidth,nHeight,saveFileName);
}