#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"


class Pedigree;

class LinkageProgram : public Program {
    
    bool run_pedigree(Pedigree& p);

 public :
    LinkageProgram(char* ped, char* map, char* dat, bool verbose) : 
        Program(ped, map, dat, verbose) {}
    
	~LinkageProgram() {}

    bool run();
};

#endif

