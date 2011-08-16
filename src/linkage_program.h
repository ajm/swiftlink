#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "program.h"


class Pedigree;
class Peeler;

class LinkageProgram : public Program {
    
    string output_filename;
    
    Peeler* run_pedigree(Pedigree& p);
    void free_peelers(vector<Peeler*>& p);

 public :
    LinkageProgram(char* ped, char* map, char* dat, char* outputfile, bool verbose) : 
        Program(ped, map, dat, verbose), 
        output_filename(outputfile) {}
    
	~LinkageProgram() {}
    
    bool run();
};

#endif

