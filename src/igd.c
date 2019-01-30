//=====================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//
//01/24/2019: Change definition of a region/interval to be half-open
//	Implement AIList for in-bin search: an option
//	Test with new data
//-------------------------------------------------------------------------------------

#include "igd_base.h"
#include "igd_create.h"
#include "igd_search.h"
#define PROGRAM_NAME  "igd"
#define MAJOR_VERSION "0"
#define MINOR_VERSION "1"
#define REVISION_VERSION "1"
#define BUILD_VERSION "0"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION

char *fileBase = "_b14_";         
uint32_t nmax[] = {15940, 15580, 12760, 12240, 11670, 10990, 10260, 9370, 8860, 8610, 8710, 
        8580, 7300, 6840, 6510, 5830, 5370, 5160, 3820, 4150, 2980, 3240, 9880, 3510};
char *folder[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
     	 "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
uint32_t gstart[25] = {0, 15940, 31520, 44280, 56520, 68190, 79180, 89440, 98810, 107670, 
        116280, 124990, 133570, 140870, 147710, 154220, 160050, 165420, 170580, 174400, 
        178550, 181530, 184770, 194650, 198160}; 
uint32_t nTiles = 198160;
uint32_t nbp = 16384, bgz_buf = 1024;
uint64_t maxCount = 268435456;//*16 Bytes;//536870912;		
uint32_t *g2ichr;

int igd_help(int argc, char **argv, int exit_code);

int main(int argc, char **argv)
{
    if (argc < 4 || argc > 8) return igd_help(argc, argv, 0);
    char *cmd = argv[1];

    if (strcmp(cmd,"create") == 0 && argc >= 5){
        return igd_create(argc, argv);
    }
    else if (strcmp(cmd,"search") == 0 && argc >= 4){
        return igd_search(argc, argv);
    }
    else {
        fprintf(stderr, "Unknown command\n");
        return igd_help(argc, argv, EX_USAGE);
    }
}

int igd_help(int argc, char **argv, int exit_code)
{
    fprintf(stderr,
"%s, v%s\n" 
"usage:   %s <command> [options]\n" 
"         create    Create an igd database\n"
"         search    Search an igd database\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}
