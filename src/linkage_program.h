#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"


class Pedigree;
class Peeler;

class LinkageProgram : public Program {
    
    string output_filename;
    unsigned int iterations;
    
    Peeler* run_pedigree(Pedigree& p);
    void free_peelers(vector<Peeler*>& p);

 public :
    LinkageProgram(char* ped, char* map, char* dat, char* outputfile, unsigned int iterations, bool verbose) : 
        Program(ped, map, dat, verbose), 
        output_filename(outputfile),
        iterations(iterations) {}
    
	~LinkageProgram() {}
    
    bool run();
};

#endif

