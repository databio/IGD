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

int igd_help(int argc, char **argv, int exit_code);
void *hc;				//extern from igd_base.h
iGD_t *IGD;
gdata_t *gData = NULL;
int32_t preIdx, preChr;
FILE *fP;

int main(int argc, char **argv)
{
    if (argc < 2) return igd_help(argc, argv, 0);
    char *cmd = argv[1];

    if (strcmp(cmd, "create") == 0){
        return igd_create(argc, argv);
    }
    else if (strcmp(cmd, "search") == 0){ 	
        return igd_search(argc, argv);
        if(gData!=NULL)free(gData);
        if(fP!=NULL)fclose(fP);
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
