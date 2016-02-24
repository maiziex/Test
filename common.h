#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>

#define VCF_MANDATORY 8
#define MAX_CHROMOSOME_LENGTH 3e6


FILE* fileOpenR(char* filename);
FILE* fileOpenRB(char* filename);
FILE* fileOpenW(char* filename);
FILE* fileOpenWB(char* filename);
FILE* fileOpenA(char* filename);
void fileClose(FILE* file);
int fileExists(const char *filename);
#endif
