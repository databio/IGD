//=====================================================================================
//Read igd region data and query data, and then find all overlaps
//by Jianglin Feng  05/12/2018
//
//database intervals sorted by _start: 8/12/2019
//Reorg--simplify: 11/06/19
//-------------------------------------------------------------------------------------
#include "igd_search.h"
//-------------------------------------------------------------------------------------
int search_help(int exit_code)
{
    fprintf(stderr,
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

void get_overlaps(iGD_t *iGD, char *chrm, int32_t qs, int32_t qe, int64_t *hits)
{
	int ichr = get_id(iGD, chrm);
	if(ichr<0)
		return;
	int i, j, n1 = qs/iGD->nbp, n2 = (qe-1)/iGD->nbp;	//define boundary!
	int32_t tE, tS, tL, tR, tM, tmpi, tmpi1, mlen, mTile = iGD->nTile[ichr]-1;
	if(n1>mTile)
		return;
	n2 = MIN(n2, mTile);
	tmpi = iGD->nCnt[ichr][n1];
	tmpi1 = tmpi-1;
	if(tmpi>0){
		if(n1!=iGD->preIdx || ichr!=iGD->preChr){
			fseek(iGD->fP, iGD->tIdx[ichr][n1], SEEK_SET);
			free(iGD->gData);
			iGD->gData = malloc(tmpi*sizeof(gdata_t));
			fread(iGD->gData, sizeof(gdata_t)*tmpi, 1, iGD->fP);
			iGD->preIdx = n1;
			iGD->preChr = ichr;
		}
		if(qe>iGD->gData[0].start){						//sorted by start
			//find the 1st rs < qe
			tL = 0, tR=tmpi1;
			while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
				tM = (tL+tR)/2;
				if(iGD->gData[tM].start < qe)	//right side:
				    tL = tM;
				else
				    tR = tM;				//left side
			}
			if(iGD->gData[tR].start<qe)tL = tR;
			//-----------------------------------------
			for(i=tL; i>=0; i--){
				if(iGD->gData[i].end>qs){
					hits[iGD->gData[i].idx]++;
				}
			}
		}
		if(n2>n1){									//n2>n1
			int32_t bd = iGD->nbp*(n1+1);			//only keep the first
			for(j=n1+1; j<=n2; j++){				//n2 inclusive!!!
				tmpi = iGD->nCnt[ichr][j];
				tmpi1 = tmpi-1;
				if(tmpi>0){
					if(j!=iGD->preIdx || ichr!=iGD->preChr){
						fseek(iGD->fP, iGD->tIdx[ichr][j], SEEK_SET);
						free(iGD->gData);
						iGD->gData = malloc(tmpi*sizeof(gdata_t));
						fread(iGD->gData, sizeof(gdata_t)*tmpi, 1, iGD->fP);
						iGD->preIdx = j;
						iGD->preChr = ichr;
					}
					if(qe>iGD->gData[0].start){
						tS = 0;
						while(tS<tmpi && iGD->gData[tS].start<bd)tS++;	//qs<bd
						tL = 0, tR=tmpi1;
						while(tL<tR-1){					//result: tR=tL+1, tL.s<qe
							tM = (tL+tR)/2;
							if(iGD->gData[tM].start < qe)	//right side:
								tL = tM;
							else
								tR = tM;				//left side
						}
						if(iGD->gData[tR].start<qe)tL = tR;
						//-----------------------------------------
						for(i=tL; i>=tS; i--){
							if(iGD->gData[i].end>qs){
								hits[iGD->gData[i].idx]++;
							}
						}
					}
				}
				bd+=iGD->nbp;
			}
		}
	}
}

void getOverlaps(char **igdFile, char **qFile, int64_t *hits)
{
	iGD_t *iGD = open_iGD(*igdFile);
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	if ((fp = gzopen(*qFile, "r")) == 0)
		return;
	ks = ks_init(fp);
    char *chrm;
	int32_t st, en, nl=0;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		chrm = parse_bed(str.s, &st, &en);
		if (chrm) {
			get_overlaps(iGD, chrm, st, en, hits);
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	close_iGD(iGD);
}

