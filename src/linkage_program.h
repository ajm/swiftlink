#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"
#include "types.h"

class Pedigree;
class Peeler;

class LinkageProgram : public Program {
    
    string output_filename;
    struct mcmc_options options;
    bool use_gpu;
    
    double* run_pedigree(Pedigree& p);
    void free_lodscores(vector<double*>& p);

 public :
    LinkageProgram(char* ped, char* map, char* dat, char* outputfile, struct mcmc_options options, bool verbose) : 
        Program(ped, map, dat, verbose),
        output_filename(outputfile),
        options(options),
        use_gpu(false) {}
    
	~LinkageProgram() {}
    
    bool run();
    
    void set_usegpu(bool _use_gpu) {
        use_gpu = _use_gpu;
    }
};

#endif

