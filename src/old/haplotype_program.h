#ifndef LKG_HAPLOTYPEPROGRAM_H_
#define LKG_HAPLOTYPEPROGRAM_H_

#include "program.h"


class HaplotypeProgram : public Program {
    
 public :
    HaplotypeProgram(char* ped, char* map, char* dat, bool verbose)
		: Program(ped, map, dat, verbose) {}
    
	~HaplotypeProgram() {}

    bool run();
    bool run_pedigree_sa(Pedigree& p);
};

#endif

