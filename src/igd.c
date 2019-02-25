//=====================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//
//01/24/2019: Change definition of a region/interval to be half-open, 0-based
//	Implement AIList for in-bin search: an option
//	Test with new data
//-------------------------------------------------------------------------------------
#include "igd_base.h"
#include "igd_create.h"
#include "igd_search.h"

char *fileBase = "_b14_";         
uint32_t nmax[] = {16040, 15680, 12860, 12340, 11770, 11090, 10360, 9470, 
        8960, 8710, 8810, 8680, 7400, 6940, 6610, 5930, 
        5470, 5260, 3920, 4250, 3080, 3340, 9980, 3710,
        2, 302, 300, 296, 290, 284, 284, 284, 106, 38, 
        36, 14, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 10, 10, 
        10, 10, 10, 10, 8, 8, 8, 6, 6, 4, 
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 2, 2, 2, 2, 2};
char *folder[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
     	 "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", 
     	 "chrM", "chr6_ssto_hap7", "chr6_mcf_hap5", "chr6_cox_hap2", "chr6_mann_hap4", "chr6_apd_hap1",
     	 "chr6_qbl_hap6", "chr6_dbb_hap3", "chr17_ctg5_hap1", "chr4_ctg9_hap1", "chr1_gl000192_random",
 	    "chrUn_gl000225", "chr4_gl000194_random", "chr4_gl000193_random", "chr9_gl000200_random", 
 	    "chrUn_gl000222", "chrUn_gl000212", "chr7_gl000195_random", "chrUn_gl000223", "chrUn_gl000224",
 	    "chrUn_gl000219", "chr17_gl000205_random", "chrUn_gl000215", "chrUn_gl000216", "chrUn_gl000217",
 	    "chr9_gl000199_random", "chrUn_gl000211", "chrUn_gl000213", "chrUn_gl000220", "chrUn_gl000218",
 	    "chr19_gl000209_random", "chrUn_gl000221", "chrUn_gl000214", "chrUn_gl000228", "chrUn_gl000227",
 	    "chr1_gl000191_random", "chr19_gl000208_random", "chr9_gl000198_random", "chr17_gl000204_random", 
 	    "chrUn_gl000233", "chrUn_gl000237", "chrUn_gl000230", "chrUn_gl000242", "chrUn_gl000243", 
 	    "chrUn_gl000241", "chrUn_gl000236", "chrUn_gl000240", "chr17_gl000206_random", "chrUn_gl000232",
 	    "chrUn_gl000234", "chr11_gl000202_random", "chrUn_gl000238", "chrUn_gl000244", "chrUn_gl000248", 
 	    "chr8_gl000196_random", "chrUn_gl000249", "chrUn_gl000246",  "chr17_gl000203_random", 
 	    "chr8_gl000197_random", "chrUn_gl000245", "chrUn_gl000247", "chr9_gl000201_random", "chrUn_gl000235",
 	    "chrUn_gl000239", "chr21_gl000210_random", "chrUn_gl000231","chrUn_gl000229", "chrUn_gl000226", 
 	    "chr18_gl000207_random"};
uint32_t gstart[] = {0, 16040, 31720, 44580, 56920, 68690, 79780, 90140, 99610, 108570,
        117280, 126090, 134770, 142170, 149110, 155720, 161650, 167120, 172380, 176300, 
        180550, 183630, 186970, 196950, 200660, 200662, 200964, 201264, 201560, 201850,
        202134, 202418, 202702, 202808, 202846, 202882, 202896, 202908, 202920, 202932,
        202944, 202956, 202968, 202980, 202992, 203004, 203016, 203028, 203040, 203052,
        203064, 203076, 203088, 203098, 203108, 203118, 203128, 203138, 203148, 203156,
        203164, 203172, 203178, 203184, 203188, 203192, 203196, 203200, 203204, 203208,
        203212, 203216, 203220, 203224, 203228, 203232, 203236, 203240, 203244, 203248,
        203252, 203256, 203260, 203264, 203268, 203272, 203276, 203280, 203284, 203286,
        203288, 203290, 203292, 203294}; 
uint32_t nTiles = 203400;

/*/
uint32_t nmax[] = {15940, 15580, 12760, 12240, 11670, 10990, 10260, 9370, 8860, 8610, 8710, 
        8580, 7300, 6840, 6510, 5830, 5370, 5160, 3820, 4150, 2980, 3240, 9880, 3510};
char *folder[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
     	 "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
uint32_t gstart[25] = {0, 15940, 31520, 44280, 56520, 68190, 79180, 89440, 98810, 107670, 
        116280, 124990, 133570, 140870, 147710, 154220, 160050, 165420, 170580, 174400, 
        178550, 181530, 184770, 194650, 198160}; 
uint32_t nTiles = 198160;
*/

uint32_t nbp = 16384, bgz_buf = 1024;
uint64_t maxCount = 268435456;//134217728;//*16 Bytes;//536870912;		
uint32_t *g2ichr;

int igd_help(int argc, char **argv, int exit_code);

int main(int argc, char **argv)
{
    if (argc < 2) return igd_help(argc, argv, 0);
    char *cmd = argv[1];

    if (strcmp(cmd,"create") == 0){
        return igd_create(argc, argv);
    }
    else if (strcmp(cmd,"search") == 0){
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
