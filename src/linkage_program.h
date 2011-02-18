#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"


class LinkageProgram : public Program {
    
 public :
    LinkageProgram(char* ped, char* map, char* dat)
		: Program(ped, map, dat) {}
    
	~LinkageProgram() {}

    bool run();

};

#endif

