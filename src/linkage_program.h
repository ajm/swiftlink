#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"
#include "types.h"

class Pedigree;
class Peeler;
class LODscores;

class LinkageProgram : public Program {
    
    LODscores* run_pedigree(Pedigree& p);
    LODscores* run_pedigree_average(Pedigree& p, int repeats);

 public :
    LinkageProgram(char* ped, char* map, char* dat, char* outputfile, struct mcmc_options options) : 
        Program(ped, map, dat, outputfile, options) {}
    
	~LinkageProgram() {}
    
    bool run();
};

#endif

