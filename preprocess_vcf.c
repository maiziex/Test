#include "string.h"
#include "stdlib.h"
#include "limits.h"
#include "math.h"
#include "common.h"

#define DIST_PAR 5e-5

typedef struct
{
    int64_t candidate_num;
    int64_t max_candidate_num;
    int64_t* candidate_loci;
    int* flag;
    int** GT;
    float* AF;
    int* snporindel_len;
    
} block_t;

typedef struct
{
    int64_t sample_num;
    int sample_name_max_length;
} pars_t;

typedef struct
{
    int64_t candidate_num;
    int64_t max_candidate_num;
    int64_t* candidate_loci;
    float* knn;
    float* AF;
    int* snporindel_len;
    int* flag;
    int** GT;
} indel_knn_t;



block_t* init_block(int64_t max_candidate_num, block_t* block, pars_t* pars)
{
    block->candidate_num=0;
    block->max_candidate_num = max_candidate_num;
    block->candidate_loci = (int64_t*) malloc(block->max_candidate_num*sizeof(int64_t));
    block->flag = (int*) malloc(block->max_candidate_num*sizeof(int));
    block->snporindel_len = (int*) malloc(block->max_candidate_num*sizeof(int));
    block->AF = (float*) malloc(block->max_candidate_num*sizeof(float));
    block->GT = (int**) malloc(block->max_candidate_num*sizeof(int*));
    int i;
    for(i=0;i<block->max_candidate_num;i++)
    {
        block->GT[i]=(int*) malloc(pars->sample_num*sizeof(int));
    }
    return block;
}

indel_knn_t* init_indel_knn(int64_t max_candidate_num,indel_knn_t* indel_knn,pars_t* pars)
{
    indel_knn->candidate_num = 0;
    indel_knn->max_candidate_num = max_candidate_num;
    indel_knn->flag = (int*) malloc(indel_knn->max_candidate_num*sizeof(int));
    indel_knn->candidate_loci= (int64_t*) malloc(indel_knn->max_candidate_num*sizeof(int64_t));
    indel_knn->knn= (float*) malloc(indel_knn->max_candidate_num*sizeof(float));
    indel_knn->AF= (float*) malloc(indel_knn->max_candidate_num*sizeof(float));
    indel_knn->snporindel_len= (int*) malloc(indel_knn->max_candidate_num*sizeof(int));
    indel_knn->GT = (int**) malloc(indel_knn->max_candidate_num*sizeof(int*));
    int i;
    for(i=0;i<indel_knn->max_candidate_num;i++)
    {
        indel_knn->GT[i]=(int*) malloc(pars->sample_num*sizeof(int));
    }

    return indel_knn;
}

void free_block(block_t* block)
{
    int i;
    for(i=0;i<block->max_candidate_num;i++)
    {
        free(block->GT[i]);
    }
    free(block->candidate_loci);
    free(block->flag);
    free(block->snporindel_len);
    free(block->AF);
    free(block->GT);
    free(block);
}

void free_indel_knn(indel_knn_t* indel_knn)
{
    int i;
    for(i=0;i<indel_knn->max_candidate_num;i++)
    {
        free(indel_knn->GT[i]);
    }

    free(indel_knn->candidate_loci);
    free(indel_knn->knn);
    free(indel_knn->AF);
    free(indel_knn->flag);
    free(indel_knn->snporindel_len);
    free(indel_knn->GT);
    free(indel_knn);
}


void check_recalloc_indel_knn(indel_knn_t* indel_knn,pars_t* pars,int64_t locus)
{

    indel_knn->candidate_loci[indel_knn->candidate_num] = locus;
    indel_knn->candidate_num++;
    if(indel_knn->candidate_num>= indel_knn->max_candidate_num)
    {
        indel_knn->candidate_loci = (int64_t*) realloc(indel_knn->candidate_loci, indel_knn->max_candidate_num*2*sizeof(int64_t));
        indel_knn->flag = (int*) realloc(indel_knn->flag, indel_knn->max_candidate_num*2*sizeof(int));
        indel_knn->knn = (float*) realloc(indel_knn->knn, indel_knn->max_candidate_num*2*sizeof(float));
        indel_knn->AF = (float*) realloc(indel_knn->AF, indel_knn->max_candidate_num*2*sizeof(float));
        indel_knn->snporindel_len = (int*) realloc(indel_knn->snporindel_len, indel_knn->max_candidate_num*2*sizeof(int));
        indel_knn->GT = (int**) realloc(indel_knn->GT, indel_knn->max_candidate_num*2*sizeof(int*));
        int64_t i;
        for(i = indel_knn->max_candidate_num; i < 2 * indel_knn->max_candidate_num; i++)
        {
             indel_knn->GT[i] = (int*) malloc(pars->sample_num * sizeof(int));
        }
        
        indel_knn->max_candidate_num *=2;

    }
}




void recalloc_gl_block(block_t* block, pars_t* pars)
{
    block->candidate_loci = (int64_t*) realloc(block->candidate_loci, block->max_candidate_num*2*sizeof(int64_t));
    block->flag = (int*) realloc(block->flag, block->max_candidate_num*2*sizeof(int));
    block->snporindel_len = (int*) realloc(block->snporindel_len, block->max_candidate_num*2*sizeof(int));
    block->AF = (float *) realloc(block->AF, block->max_candidate_num*2*sizeof(float));
    block->GT = (int**) realloc(block->GT, block->max_candidate_num*2*sizeof(int*));
    int64_t i;
    for(i = block->max_candidate_num; i < 2 * block->max_candidate_num; i++)
    {
        block->GT[i] = (int*) malloc(pars->sample_num * sizeof(int));
    }
    
    block->max_candidate_num *= 2;
    //printf("block size increase to %d\n", block->max_candidate_num);
}

void add_candidate_into_gl_block(block_t* block, pars_t* pars, int64_t locus)
{
    block->candidate_loci[block->candidate_num] = locus;
    block->candidate_num++;
    if(block->candidate_num >= block->max_candidate_num)
    {
        recalloc_gl_block(block, pars);
    }
}



void init_query_samples(char* input_filename, pars_t* pars)
{
    FILE* input = fileOpenR(input_filename);
    
    int field_index = 0, sample_name_length;
    pars->sample_name_max_length = 0;
    char c, c_prev;
    
    while((c = (char) getc(input)) != EOF)
    {
        if(c == '#')
        {
            c = (char) getc(input);
            if(c == '#')
            {
                while(c != '\n') c = (char) getc(input);
            }
            else
            {
                while(c != '\n')
                {
                    c_prev = c;
                    c = (char) getc(input);
                    if(c == '\n') break;
                    if((c_prev == ' ' || c_prev == '\t') && (c != ' ' && c != '\t'))
                    {
                        field_index++;
                        if(field_index > VCF_MANDATORY)
                        {
                            sample_name_length = 1;
                        }
                    }
                    else if((c_prev != ' ' || c_prev != '\t') && (c != ' ' && c != '\t'))
                    {
                        if(field_index > VCF_MANDATORY)
                        {
                            sample_name_length++;
                        }
                    }
                    else
                    {
                        if(field_index > VCF_MANDATORY && sample_name_length > pars->sample_name_max_length)
                        {
                            pars->sample_name_max_length = sample_name_length;
                        }
                    }
                }
                break;
            }
        }
    }
    
    pars->sample_num = field_index - VCF_MANDATORY;
    fileClose(input);
}




float* descend_sort(float* saved_buffer,int64_t candidate_num)
{
    int i,j,a;
    for(i=0;i<candidate_num;++i)
    {
        for(j=i+1;j<candidate_num;++j)
        {
            if(saved_buffer[i]<saved_buffer[j])
            {
                a=saved_buffer[i];
                saved_buffer[i]=saved_buffer[j];
                saved_buffer[j]=a;
            }
        }
    }
    return saved_buffer;
}




int cmpfunc(const void* a, const void*b)
{
    return (*(float*)b-*(float*)a);
}




void find_KNN2(block_t* block, pars_t* pars, indel_knn_t* indel_knn)
{
    int i,j,k,m,l,K_number;
    K_number=4;
    int union_counter;
    int intersection_counter;
    float SI;
    float knnSI;
    float* knnbucket = calloc(K_number+1,sizeof(float));
    double dist;
    int tmp1,tmp2;
   
    printf("%i\n",indel_knn->candidate_num);
    for(i = 0; i < indel_knn->candidate_num; i++)
    {
        printf("%i\n",i);
            // initialize knnbucket to zero
            knnbucket[0]=0.0;
            knnbucket[1]=0.0;
            knnbucket[2]=0.0;
            knnbucket[3]=0.0;
            knnbucket[4]=0.0;

            knnSI = 0.0;
            for(j=0; j< block->candidate_num;j++)
            {
                    union_counter=0;
                    intersection_counter=0;
                    for(k=0; k< pars->sample_num;k++)
                    {
                        tmp1=(block->GT[i][k]+1)>>1;
                        tmp2=(block->GT[j][k]+1)>>1;
                        union_counter+=(tmp1|tmp2);
                        intersection_counter+=(tmp1&tmp2);
                        /*
                        if ((indel_knn->GT[i][k]==1 || indel_knn->GT[i][k]==2) ||  (block->GT[j][k]==1 || block->GT[j][k]==2))
                        {
                            union_counter ++;
                        }
                        
                        if ((indel_knn->GT[i][k]==1 || indel_knn->GT[i][k]==2) &&  (block->GT[j][k]==1 || block->GT[j][k]==2))
                        {
                            intersection_counter ++;
                        }
                         */
                    }
                    if(intersection_counter==0)
                    {
                             continue;
                    }
                    dist = exp(-DIST_PAR*abs(indel_knn->candidate_loci[i]-block->candidate_loci[j]));
                    SI =(float) intersection_counter*(1.0)/union_counter;
                    knnbucket[4] = SI;
                   // qsort(knnbucket,5,sizeof(float),cmpfunc);
                    descend_sort(knnbucket,K_number+1);

                
            }
            
            for(l=0; l< K_number;l++)
            {
                knnSI =  knnSI + knnbucket[l];
            }
            
            indel_knn->knn[i] = knnSI/K_number;
        
    }
    free(knnbucket);
    
}



void print_indel_info(char* output_prefix, indel_knn_t* indel_knn)
{
    char* output_filename = (char*) malloc(strlen(output_prefix)+20);
    sprintf(output_filename,"%s.final",output_prefix);
    FILE* output;
    output = fileOpenW(output_filename);
    int i;
    for(i=0;i<indel_knn->candidate_num;i++)
    {
       // printf("here:%i",i);
        fprintf(output,"%i\t%i\t%i\t%f\t%f\t\n",indel_knn->candidate_loci[i],indel_knn->flag[i],indel_knn->snporindel_len[i],indel_knn->AF[i],indel_knn->knn[i]);
        printf("%i\t%i\t%i\t%f\t%f\t\n",indel_knn->candidate_loci[i],indel_knn->flag[i],indel_knn->snporindel_len[i],indel_knn->AF[i],indel_knn->knn[i]);
    }
    
 free(output_filename);
    fileClose(output);
}




void read_vcf_data(char* input_filename,char* output_prefix, pars_t* pars)
{
    init_query_samples(input_filename,pars);
    char c,c_prev,c_next;
    block_t* block =  (block_t*) calloc(1,sizeof(block_t));
    indel_knn_t* indel_knn =  (indel_knn_t*) calloc(1,sizeof(indel_knn_t));
    init_block(10000,block,pars);
    init_indel_knn(10000,indel_knn,pars);
    
    
    FILE* input = fileOpenR(input_filename);
    while( (c= (char) getc(input)) != EOF)
    {    if(c == '#')
    {
        c = (char) getc(input);
        if(c == '#')
        {
            while(c!='\n') c = (char) getc(input);
        }
        else
        {
            while(c!='\n') c= (char) getc(input);
            break;
        }
    }
    }
    
    int i, j;
    int field_length,ref_length,alt_length;
    int64_t locus;
    float af;
    int num_of_semicolon_before_AF,num_of_colon_before_GT,freeze;
    char* chr = calloc(MAX_CHROMOSOME_LENGTH, sizeof(char));
    char* ref = calloc(MAX_CHROMOSOME_LENGTH, sizeof(char));
    char* alt = calloc(MAX_CHROMOSOME_LENGTH, sizeof(char));
    c = (char) getc(input);
    
    while(c != EOF)
    {
        memset(chr, 0, MAX_CHROMOSOME_LENGTH * sizeof(char));
        field_length = 0;
        while(c != '\t')
        {
            chr[field_length] = c;
            c = (char) getc(input);
            field_length++;
        }
        
        locus = 0;
        while((c = (char) getc(input)) != '\t')
        {
            locus = locus * 10 + c - '0';
        }
        
        printf("%s:%i\n", chr, locus);
        /*if(locus==260315)
         {
         printf("hereerror\n");
         }*/
        while((c = (char) getc(input)) != '\t');
        
        memset(ref, 0, MAX_CHROMOSOME_LENGTH * sizeof(char));
        field_length = 0;
        while((c = (char) getc(input)) != '\t')
        {
            ref[field_length] = c;
            field_length++;
        }
        ref_length = field_length;
        
        
        memset(alt, 0, MAX_CHROMOSOME_LENGTH * sizeof(char));
        field_length = 0;
        while((c = (char) getc(input)) != '\t')
        {
            alt[field_length] = c;
            field_length++;
        }
        alt_length = field_length;
        
        if(ref_length == 1 && alt_length ==1)
        {
            block->flag[block->candidate_num] = 0; //snp
        }
        else
        {
            if(ref_length> alt_length)
            {
                  indel_knn->flag[indel_knn->candidate_num] = 1; //delete
                  indel_knn->snporindel_len[indel_knn->candidate_num] = abs(alt_length-ref_length);
            }
            else
            {
                indel_knn->flag[indel_knn->candidate_num]=2; //insert
                indel_knn->snporindel_len[indel_knn->candidate_num] = abs(alt_length-ref_length);
            }
        }
        
        
        
        for(i = 0; i < 2; i++)
        {
            while((c = (char) getc(input)) != '\t');
        }
        
        
        freeze=0;
        num_of_semicolon_before_AF=0;
        c = (char) getc(input);
        char* fbuffer = calloc(10,sizeof(char));
        while(c != '\t')
        {
            c_prev = c;
            c = (char) getc(input);
            if(c_prev == ';' && !freeze) num_of_semicolon_before_AF++;
            if(c_prev == ';' && c == 'A' && ((c = (char) getc(input))=='F') && ((c = (char) getc(input))=='='))
                freeze = 1;
            
            i = 0;
            if(freeze==1)
            {
                
                while((c = (char) getc(input)) != ';' && (c != '\t') && (c!= '\n'))
                {
                    fbuffer[i] = c;
                    //   printf("%c",c);
                    i++;
                }
                // printf("\n");
                af = atof(fbuffer);
                if(ref_length == 1 && alt_length ==1)  //snp)
                {
                      block->AF[block->candidate_num]= af;
                }
                else
                {
                     indel_knn->AF[indel_knn->candidate_num]=af;
                }
                freeze=0;
                
            }
        }
        free(fbuffer);
        
        ///////////////////////////////////////////////
        // read field with GT:
        freeze=0;
        num_of_colon_before_GT = 0;
        c = (char) getc(input);
        while(c != '\t')
        {
            c_prev = c;
            c = (char) getc(input);
            if(c_prev == ':'&& freeze==0) num_of_colon_before_GT++;
            if(c_prev == 'G' && c == 'T')
            {
                freeze = 1;
            }
            
        }
        
        // read field after GT
        for(i = 0; i < pars->sample_num; i++)
        {
            for(j = 0; j < num_of_colon_before_GT; j++)
            {
                while(c != ':') c = (char) getc(input);
                c = (char) getc(input);
            }
            
            c = (char) getc(input);
            while(c != ':' && c != '\t' && c != '\n')
            {
                c_prev = c;
                c = (char) getc(input);
                c_next = (char) getc(input);
       
                if(ref_length==1 && alt_length==1) //snp
                {
                       block->GT[block->candidate_num][i] = (c_prev-'0')+(c_next-'0');
                }
                else
                { 
                       indel_knn->GT[indel_knn->candidate_num][i]= (c_prev-'0')+(c_next-'0');
                }
                c = (char) getc(input);
            }
            while((c = (char) getc(input)) != '\t')
            {
                if (c == EOF)
                {
                    printf("End of file\n");
                    break;
                }
            }
            
            ///////////////////////////////////////////////
            
        }
        
        if(ref_length==1 && alt_length==1) //snp
       {
             add_candidate_into_gl_block(block,pars,locus);
       }
       else
       {
              check_recalloc_indel_knn(indel_knn,pars,locus);
       }
    }
    
    find_KNN2(block,pars,indel_knn);
   //  find_KNN(block,pars,indel_knn);
   
    printf("%i",indel_knn->candidate_num);
    printf("\n");
    print_indel_info(output_prefix, indel_knn); 
    
    free_block(block);
    free_indel_knn(indel_knn);
    free(chr);
    free(ref);
    free(alt);
}


int main(int argc, char* argv[])
{
    pars_t* pars =  (pars_t*) calloc(1,sizeof(pars_t));
    // set_default_pars(pars);
    read_vcf_data(argv[1], argv[2],pars);
    printf("DONE!!!\n");
    
    free(pars); 
    
    
    
}
