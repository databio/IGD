//=====================================================================================
//Common igd struct, parameters, functions
//by Jianglin Feng  05/12/2018
//
//01/24/2019: Change definition of a region/interval to be half-open, 0-based
//	Implement AIList for in-bin search: an option
//	Test with new data
//-------------------------------------------------------------------------------------
#include "libigd/igd_base.h"
#include "libigd/igd_create.h"
#include "libigd/igd_search.h"

// A few internal function prototypes, should probably move to an igd.h file
int igd_help(int argc, char **argv, int exit_code);
int create_help(int exit_code);
int search_help(int exit_code);

/**
 * @brief Parse the command-line args (argv) if given 'search' command
 *
 * @param argc Arg count directly from main()
 * @param argv Arg vars directly from main()
 * @return A SearchParams_t* with all parameters needed for requested task.
 */
SearchParams_t *parse_search_args(int argc, char **argv);

/**
 * @brief Parse the command-line args (argv) if given 'create' command
 *
 * @param argc Arg count directly from main()
 * @param argv Arg vars directly from main()
 * @return A CreateParams_t* with all parameters needed for requested task.
 */
CreateParams_t *parse_create_args(int argc, char **argv);


// Old globals should be removed now with the C library approach
// void *hc;				//extern from igd_base.h
// IGD_t *IGD;
// gdata_t *gData = NULL;
// gdata0_t *gData0 = NULL;
// int32_t preIdx, preChr, tile_size;
// FILE *fP;

int main(int argc, char **argv)
{
    if (argc < 2) return igd_help(argc, argv, 0);
    char *cmd = argv[1];
    int result;

    if (strcmp(cmd, "create") == 0) {
        CreateParams_t* cParams = parse_create_args(argc, argv);
        result = create_IGD_from_params(cParams);
        free(cParams);
    } else if (strcmp(cmd, "search") == 0) {
        SearchParams_t* sParams = parse_search_args(argc, argv);
        if (sParams->status == FAILED) {
            printf("Status: {%d}\n", sParams->status); 
            // free(sParams);
            return EX_OK;
        }
        // create IGD_t object
        IGD_t *IGD = open_IGD(sParams->igdFileName);
        // allocate response vector
        int32_t *hits = calloc(IGD->nFiles, sizeof(int32_t));
        // perform search
        result = getOverlapsFile(IGD, sParams->queryFileName, hits);
        // result = search_IGD(IGD, sParams);

         // write results to stdout
        printf("index\ttotal_regions\toverlaps\tfile_name\n");        
        for(int i=0;i<IGD->nFiles;i++) {
            printf("%i\t%i\t%lld\t%s\n", i, IGD->finfo[i].nr, (long long)hits[i], IGD->finfo[i].fileName);
        }
        //Destruct
        free(hits);
        free(sParams);
        close_IGD(IGD);
    } else {
        fprintf(stderr, "Unknown command\n");
        return igd_help(argc, argv, EX_USAGE);;
    }
    return EX_OK;
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

int create_help(int exit_code)
{
    printf(
"%s, v%s\n"
"usage:   %s create <input dir> <output dir> <output igd name> [options] \n"
"             -b  <Tile size in power of 2 (default 14)> \n"
"             -c  < .BED column as value >=4 (default 4) \n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

int search_help(int exit_code)
{
    printf(
"%s, v%s\n"
"usage:   %s search <igd database file> [options]\n"
"         options:\n"
"             -q <query file>\n"
"             -r <a region: chrN start end>\n"
"             -v <signal value 0-1000>\n"
"             -o <output file Name>\n"
"             -c display all intersects\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

SearchParams_t *parse_search_args(int argc, char **argv) {
    SearchParams_t *sParams = SearchParams_init();  // Allocate memory for search settings struct
    sParams->status = FAILED;  // pre-initialize to FAILED

    if(argc<4) {
        search_help(EX_OK);
        return sParams;
    }

    char* igdName = argv[2];
    int32_t v = 0, qs=1, qe=2;
    int i, i1, j, checking=0, mode=-1, mt=0, ext=0, xlen=0,  mv=0, mx=0, ichr, k;
    char out[64]="";  
    char *chrm;    
    char *qfName = "";  
    uint64_t ols;
    // Make sure given database has correct extension.
    char *ftype = igdName + strlen(igdName) - 4;
    if(strcmp(".igd", ftype)!=0){
        printf("%s is not an igd database", igdName);
        return sParams;
    }    
    FILE* fi = fopen(igdName, "rb");
    if(!fi){
        printf("%s does not exist", igdName);
        return sParams;
    }
    fclose(fi);

    for(int i=3; i<argc; i++){
        if(strcmp(argv[i], "-q")==0){
            if(i+1<argc){
                qfName = argv[i+1];
                mode = 1;
            }
            else{
                printf("No query file.\n");
                return sParams;
            }   
        }
        else if(strcmp(argv[i], "-r")==0){
            if(i+3<argc){
                mode = 2;
                chrm = argv[i+1];
                qs = atoi(argv[i+2]);
                qe = atoi(argv[i+3]);
            }        
        }
        else if(strcmp(argv[i], "-v")==0){//>=v
            mv = 1;
            if(i+1<argc)
                v = atoi(argv[i+1]);
        }
        else if(strcmp(argv[i], "-m")==0){
            mode = 0;
        } 
        else if(strcmp(argv[i], "-s")==0){
            mode = 3;   //seqpare
        }                  
        else if(strcmp(argv[i], "-o")==0){
            if(i+1<argc)
                strcpy(out, argv[i+1]);
        }   
        else if(strcmp(argv[i], "-c")==0){
            checking = 1;
        }                              
    }

    strcpy(sParams->igdFileName, igdName);
    strcpy(sParams->queryFileName, qfName);
    sParams->checking = checking;
    sParams->datamode = mode;
    sParams->status = SUCCESS;
    return sParams;
}  

CreateParams_t * parse_create_args(int argc, char **argv) {
    if (argc < 5) 
        return create_help(EX_OK);       

    CreateParams_t *cParams = CreateParams_init();

    int32_t i, j, n;
    // char ipath[1024];
    // char opath[1024];
    // char *dbname = argv[4]; 
    char *s1, *s2;
    strcpy(cParams->inputPath, argv[2]);
    strcpy(cParams->outputPath, argv[3]);
    strcpy(cParams->igdName, argv[4]);
    cParams->datamode = 1;
    cParams->filetype = 0;
    int dtype = 1, ftype = 0;   //file type 1: list of bed file
    cParams->tile_size = 16384;                              //2^14 = 16384; arg -b 14 (default)

    // process command-line arguments
    for(i=5; i<argc; i++){
        if(strcmp(argv[i], "-s")==0 && i+1<argc)    //data structure type
            cParams->datamode = atoi(argv[i+1]); 
        if(strcmp(argv[i], "-b")==0 && i+1<argc){
            n = atoi(argv[i+1]);
            if(n>10 && n<20)
                cParams->tile_size = pow(2, n);  
        } 
        if(strcmp(argv[i], "-f")==0) 
            cParams->filetype = 1;                                
    }
    
    if(cParams->outputPath[strlen(cParams->outputPath)-1]!='/'){
        strcat(cParams->outputPath, "/");
    }          
          
    if(cParams->filetype==0 && cParams->datamode!=2){                            
        if(cParams->inputPath[strlen(cParams->inputPath)-1]=='/'){
            strcat(cParams->inputPath, "*");
        }
        else if(cParams->inputPath[strlen(cParams->inputPath)-1]!='*'){
            strcat(cParams->inputPath, "/*");
        }
    }  
    return cParams;
}

