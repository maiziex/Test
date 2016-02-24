
#include <stdlib.h>
#include "common.h"



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


