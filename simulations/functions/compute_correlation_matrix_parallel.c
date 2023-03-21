#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int computeFileSize(char *filename, int* saveWidth, int* saveHeight)
{
    FILE* stream = fopen(filename, "r");
    if(!stream){
        printf("ERROR: %s",filename);
    }
    char line[10000];
    int nWidth = 0;
    int nHeight = 0;
    char *token;
    while (fgets(line, 10000, stream))
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
    char line[10000];
    int i = 0;
    int j = 0;
    while (fgets(line, 10000, stream))
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
    fprintf(errorFile,"%d,%d\n",sizeof(int) * nWidth*nHeight,nWidth*i+j);
    fprintf(errorFile,"%d,%d\n",5,10);
    fclose(stream);
    fclose(errorFile);

    return data;
}

int computePairwiseCorrelation(int *data,int nWidth, int nHeight, char*saveFileName)
{
    FILE* outputFile = fopen(saveFileName,"w");

    int *sparseCorrelationMatrix;
    sparseCorrelationMatrix = malloc((sizeof(int)*2+sizeof(float)) * nHeight*nHeight);


    int m = ceil(log2f(nHeight+1)-1);
    int arrayLen[m];
    for(int j = 0; j<m; j++){
        arrayLen[j] = nHeight;
        for(int i = 0; i<j; i++){
            arrayLen[j] *= 2;
        }
    }

    #pragma omp parallel
    {
        // Initialize pointers to arrays with length 2^n
        int *neuron1_local[m], *neuron2_local[m];
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
        int totalCounter;
        int arrayCounter = 0;
        int runningLength = 0;
        int totalLength = 0;

        #pragma omp for
        for (i = 0; i<nHeight; i++) {
            for (j = i+1; j<nHeight; j++) {
                i0 = i*nWidth;
                overlapCounter = 0;
                totalCounter = 0;
                while(data[i0]!=0) {
                    j0 = j*nWidth;
                    while(data[j0]!=0) {
                        if(abs(data[i0]-data[j0])<10) {
                            overlapCounter++;
                        }
                        j0++;
                    }
                    i0++;
                }
                if(overlapCounter>1){
                    totalCounter = (i0-i*nWidth)*(j0-j*nWidth);
                    R = (float)overlapCounter / (float)totalCounter;
                    if(R>=0.0001){
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
                    fprintf(outputFile,"%d,%d,%.4f\n",neuron1_local[j][i],neuron2_local[j][i],correlation_local[j][i]);
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