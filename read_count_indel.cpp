#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <unordered_map>
#include <string>
#include <cctype>
#include <stdexcept>
#include <stdio.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>


using namespace std;
#define MAX_POSSIBLILITY_LENGTH 50

unordered_map<string,long> dict_indel;

void toUpper(basic_string<char>& s) {
    for (basic_string<char>::iterator p = s.begin();
         p != s.end(); ++p) {
        *p = toupper(*p); // toupper is for char
    }
}



void read_pileup(string input_filename,  long long start, long long end)
{
    ifstream input;
    input.open(input_filename.c_str());
    
    int i, j,ii;
    
    string chr, content, dummy;
    int locus, copy, num;
    char ref;
    string current_indel_key;
    
    
    while(input >> chr >> locus >> ref >> copy)
    {
        if(copy == 0) continue;
        input >> content >> dummy;
        
        if(locus < start || locus >= end) continue;
        
        
        if(content.find("+")!= std::string::npos ||  content.find("-")!= std::string::npos)
        {
            //printf("%i\t%c\t%s\n", locus,ref,content.c_str());
            for(i = 0; i < content.size(); i++)
            
            if((content[i] == '.' || content[i] == ',') && content[i+1] != '+' && content[i+1] != '-')
            {
                if(ref == 'A' || ref == 'a')
                {
                    //printf("here 1\n");
                    current_indel_key = to_string(locus) + 'A' ;
                    auto it = dict_indel.find(current_indel_key);
                    if(it != dict_indel.end())
                    it->second = dict_indel[current_indel_key]+1 ;
                    else
                    dict_indel[current_indel_key] = 1;
                    
                }
                else if(ref == 'T' || ref == 't')
                {
                    //printf("here 2\n");
                    current_indel_key = to_string(locus) + 'T';
                    auto it = dict_indel.find(current_indel_key);
                    if(it != dict_indel.end())
                    it->second = dict_indel[current_indel_key]+1 ;
                    else
                    dict_indel[current_indel_key] = 1;
                    
                }
                
                else if(ref == 'C' || ref == 'c')
                {
                    //printf("here 3\n");
                    current_indel_key = to_string(locus) + 'C';
                    auto it = dict_indel.find(current_indel_key);
                    if(it != dict_indel.end())
                    it->second = dict_indel[current_indel_key]+1 ;
                    else
                    dict_indel[current_indel_key] = 1;
                    
                }
                
                else if(ref == 'G' || ref == 'g')
                {
                    //printf("here 4\n");
                    current_indel_key = to_string(locus) + 'G';
                    auto it = dict_indel.find(current_indel_key);
                    if(it != dict_indel.end())
                    it->second = dict_indel[current_indel_key]+1 ;
                    else
                    dict_indel[current_indel_key] = 1;
                    
                }
                
            }
            else if((content[i] == 'A' || content[i] == 'a') && content[i+1] != '+' && content[i+1] != '-')
            {
                //printf("here 5\n");
                current_indel_key = to_string(locus) + 'A';
                auto it = dict_indel.find(current_indel_key);
                if(it != dict_indel.end())
                it->second = dict_indel[current_indel_key]+1 ;
                else
                dict_indel[current_indel_key] = 1;
                
            }
            else if((content[i] == 'T' || content[i] == 't') && content[i+1] != '+' && content[i+1] != '-')
            {
                //printf("here 6\n");
                current_indel_key = to_string(locus) + 'T';
                auto it = dict_indel.find(current_indel_key);
                if(it != dict_indel.end())
                it->second = dict_indel[current_indel_key]+1 ;
                else
                dict_indel[current_indel_key] = 1;
                
            }
            else if((content[i] == 'C' || content[i] == 'c') && content[i+1] != '+' && content[i+1] != '-')
            {
                //printf("here 7\n");
                current_indel_key = to_string(locus) + 'C';
                auto it = dict_indel.find(current_indel_key);
                if(it != dict_indel.end())
                it->second = dict_indel[current_indel_key]+1 ;
                else
                dict_indel[current_indel_key] = 1;
                
            }
            else if((content[i] == 'G' || content[i] == 'g') && content[i+1] != '+' && content[i+1] != '-')
            {
                //printf("here 8\n");
                current_indel_key = to_string(locus) + 'G';
                auto it = dict_indel.find(current_indel_key);
                if(it != dict_indel.end())
                it->second = dict_indel[current_indel_key]+1 ;
                else
                dict_indel[current_indel_key] = 1;
                
            }
            else if(content[i] == 'N' || content[i] == 'n')
            {
            }
            else if(content[i] == '^')
            {
                i++;
            }
            else if(content[i] == '$')
            {
            }
            else if(content[i] == '*')
            {
            }
            else if(content[i] == '+' || content[i] == '-')
            {
                num = 0;
                for(j = i + 1; j < content.size(); j++)
                {
                    if(content[j] >= '0' && content[j] <= '9')
                    {
                        num = num * 10 + content[j] - '0';
                    }
                    
                    else
                    {
                        //printf("here 9\n");
                        current_indel_key ="";
                        for(ii=0;ii<num+3;ii++)
                        {   if(content[j-3+ii] == '.' || content[j-3+ii] == ',')
                            {
                                current_indel_key= current_indel_key + ref;
                            }
                            else
                            {
                                current_indel_key= current_indel_key + content[j-3+ii];
                            }
                        }
                        
                        toUpper(current_indel_key);
                        current_indel_key = to_string(locus) + current_indel_key;
                        // printf("%s\n", current_indel_key.c_str());
                        auto it = dict_indel.find(current_indel_key);
                        if(it != dict_indel.end())
                        it->second = dict_indel[current_indel_key]+1 ;
                        else
                        dict_indel[current_indel_key] = 1;
                        
                        i = j + num - 1;
                        break;
                    }
                }
            }
            else
            {
                //printf("warning: unknown character %c\n", content[i]);
            }
        }
    }
    printf("finished\n");
    
    cout<<dict_indel["43990987T-8TTTATTTA"]<<endl;
    cout<<"size is "<< dict_indel.size()<<endl;
    input.close();
}

FILE* fileOpenWB(const char* filename)
{
    FILE* file = (FILE*) fopen(filename, "wb");
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

int main(int argc, char* argv[])
{
    struct dirent *pDirent;
    struct stat filestat;
    string filepath;
    DIR *pDir;
    string input_folder = argv[1];
    string output_filename = argv[2];
    string start_locus = argv[3];
    string end_locus = argv[4];
    long long start = atoll(start_locus.c_str());
    long long end = atoll(end_locus.c_str());
    long long genome_len = end - start;
    
    
    
    pDir = opendir (argv[1]);
    if (pDir == NULL)
    {
        //printf ("Cannot open directory '%s'\n", input_folder);
        cout << "Error(" << errno << ") opening " << pDir << endl;
        return 0;
    }
    
    while ((pDirent = readdir(pDir)) != NULL)
    {
        filepath= input_folder + "/" + pDirent->d_name;
        
        if (stat( filepath.c_str(), &filestat ))
        continue;
        if (S_ISDIR( filestat.st_mode ))
        continue;
        
        
        read_pileup(filepath,start, end);
    }
    
    
    ofstream myfile("saved_indel.txt");
    
    if(myfile.is_open())
    {
        for(unordered_map<string, long, less< int > >::const_iterator iter = dict_indel.begin();
            iter != dict_indel.end(); ++iter )
        {
            // fprintf(output,"%s\t%lld\n",iter->first,iter->second);
            myfile<< iter->first << '\t' << iter->second << '\n';
        }
        
        myfile.close();
    }
    
    
    //fileClose(output);
    closedir (pDir);
    return 0;
}
