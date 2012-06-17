#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"
#include "types.h"

class Pedigree;
class Peeler;

class LinkageProgram : public Program {
    
    double* run_pedigree(Pedigree& p);
    void free_lodscores(vector<double*>& p);

 public :
    LinkageProgram(char* ped, char* map, char* dat, char* outputfile, struct mcmc_options options) : 
        Program(ped, map, dat, outputfile, options) {}
    
	~LinkageProgram() {}
    
    bool run();
};

#endif

