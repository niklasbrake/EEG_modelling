#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    int fclose(steam);
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
    FILE* errorFile = fopen("C:/Users/brake/workspace/errors.txt","w");
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
    int fclose(steam);
    int fclose(errorFile);

    return data;
}

int computePairwiseCorrelation(int *data,int nWidth, int nHeight, char*saveFileName)
{
    FILE* outputFile = fopen(saveFileName,"w");
    float R;
    int i;
    int j;
    int i0;
    int j0;
    int overlapCounter;
    int totalCounter;
    for (i = 0; i<nHeight; i++)
    {
        for (j = i+1; j<nHeight; j++)
        {
            i0 = i*nWidth;
            overlapCounter = 0;
            totalCounter = 0;
            while(data[i0]!=0)
            {
                j0 = j*nWidth;
                while(data[j0]!=0)
                {
                    if(abs(data[i0]-data[j0])<10)
                    {
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
                    fprintf(outputFile,"%d,%d,%.4f\n",i,j,R);
                }
            }
        }
    }
    int fclose(outputFile);
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