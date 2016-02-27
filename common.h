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
int chr_equal(char* chr, char* chr_prev);
void chr_copy(char* chr, char* chr_prev);
int getMaxInt(int a, int b);
int getMinInt(int a, int b);
int64_t getMinInt64(int64_t a, int64_t b);
#endif
