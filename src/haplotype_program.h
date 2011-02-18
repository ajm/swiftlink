#ifndef LKG_HAPLOTYPEPROGRAM_H_
#define LKG_HAPLOTYPEPROGRAM_H_

#include "program.h"


class HaplotypeProgram : public Program {
    
 public :
    HaplotypeProgram(char* ped, char* map, char* dat)
		: Program(ped, map, dat) {}
    
	~HaplotypeProgram();

    bool run();

};

#endif

