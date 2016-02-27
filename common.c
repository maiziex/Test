
#include <stdlib.h>
#include "common.h"

int chr_equal(char* chr, char* chr_prev)
{
        int i;
        for(i = 0; i < MAX_CHROMOSOME_LENGTH; i++)
        {
                if(chr[i] != chr_prev[i]) return 1;
        }
        return 0;
}

void chr_copy(char* chr, char* chr_prev)
{
        int i;
        for(i = 0; i < MAX_CHROMOSOME_LENGTH; i++)
        {
                chr_prev[i] = chr[i];
        }
}

int getMaxInt(int a, int b)
{
	if(a > b) return a; else return b;
}

int getMinInt(int a, int b)
{
	if(a < b) return a; else return b;
}

int64_t getMinInt64(int64_t a, int64_t b)
{
        if(a < b) return a; else return b;
}

FILE* fileOpenR(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "r");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenRB(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "rb");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenW(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "w");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenWB(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "wb");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenA(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "a");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

void fileClose(FILE* file)
{
        fclose(file);
}

int fileExists(const char *filename)
{
        FILE *file = fopen(filename, "r");
        if(file)
        {
                fclose(file);
                return 1;
        }
        return 0;
}


