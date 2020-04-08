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

SearchTask_t *parse_search_args(int argc, char **argv);
CreateTask_t *parse_create_args(int argc, char **argv);

void *hc;				//extern from igd_base.h
IGD_t *IGD;
gdata_t *gData = NULL;
// gdata0_t *gData0 = NULL;
int32_t preIdx, preChr, tile_size;
FILE *fP;

int main(int argc, char **argv)
{
    if (argc < 2) return igd_help(argc, argv, 0);
    char *cmd = argv[1];
    int result;

    if (strcmp(cmd, "create") == 0) {
        CreateTask_t* cTask = parse_create_args(argc, argv);
        result = create_IGD_from_task(cTask);
        free(cTask);
    } else if (strcmp(cmd, "search") == 0) {
        SearchTask_t* sTask = parse_search_args(argc, argv);
        if (sTask->status == FAILED) {
            printf("SearchTask status: {%d}\n", sTask->status); 
            printf("Creating task failed.\n");
            // free(sTask);
            return EX_OK;
        }
        // create IGD_t object
        IGD_t *IGD = open_IGD(sTask->igdFileName);
        // allocate response vector
        int32_t *hits = calloc(IGD->nFiles, sizeof(int32_t));
        // perform search
        result = getOverlapsFile(IGD, sTask->queryFileName, hits);
        // result = search_IGD(IGD, sTask);

         // write results to stdout
        printf("index\ttotal_regions\toverlaps\tfile_name\n");        
        for(int i=0;i<IGD->nFiles;i++) {
            printf("%i\t%i\t%lld\t%s\n", i, IGD->finfo[i].nr, (long long)hits[i], IGD->finfo[i].fileName);
        }
        //Destruct
        free(hits);   
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

SearchTask_t *parse_search_args(int argc, char **argv) {
    SearchTask_t *sTask = SearchTask_init();

    printf("Set status 1\n"); 
    sTask->status = FAILED;
    printf("Set status 2\n");
    if(argc<4) {
        search_help(EX_OK);
        return sTask;
    }
    // SearchTask_t sTask = *sTask;

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
        return sTask;
    }    
    FILE* fi = fopen(igdName, "rb");
    if(!fi){
        printf("%s does not exist", igdName);
        return sTask;
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
                return sTask;
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

    sTask->status = SUCCESS;
    strcpy(sTask->igdFileName, igdName);
    strcpy(sTask->queryFileName, qfName);
    sTask->checking = checking;
    sTask->datamode = mode;
    return sTask;
}  

CreateTask_t * parse_create_args(int argc, char **argv) {
    if (argc < 5) 
        return create_help(EX_OK);       

    CreateTask_t *cTask = CreateTask_init();

    int32_t i, j, n;
    // char ipath[1024];
    // char opath[1024];
    // char *dbname = argv[4]; 
    char *s1, *s2;
    strcpy(cTask->inputPath, argv[2]);
    strcpy(cTask->outputPath, argv[3]);
    strcpy(cTask->igdName, argv[4]);
    cTask->datamode = 1;
    cTask->filetype = 0;
    int dtype = 1, ftype = 0;   //file type 1: list of bed file
    cTask->tile_size = 16384;                              //2^14 = 16384; arg -b 14 (default)

    // process command-line arguments
    for(i=5; i<argc; i++){
        if(strcmp(argv[i], "-s")==0 && i+1<argc)    //data structure type
            cTask->datamode = atoi(argv[i+1]); 
        if(strcmp(argv[i], "-b")==0 && i+1<argc){
            n = atoi(argv[i+1]);
            if(n>10 && n<20)
                cTask->tile_size = pow(2, n);  
        } 
        if(strcmp(argv[i], "-f")==0) 
            cTask->filetype = 1;                                
    }
    
    if(cTask->outputPath[strlen(cTask->outputPath)-1]!='/'){
        strcat(cTask->outputPath, "/");
    }          
          
    if(cTask->filetype==0 && cTask->datamode!=2){                            
        if(cTask->inputPath[strlen(cTask->inputPath)-1]=='/'){
            strcat(cTask->inputPath, "*");
        }
        else if(cTask->inputPath[strlen(cTask->inputPath)-1]!='*'){
            strcat(cTask->inputPath, "/*");
        }
    }  
    return cTask;
}

