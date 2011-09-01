#ifndef LKG_TESTPROGRAM_H_
#define LKG_TESTPROGRAM_H_

#include "program.h"

class TestProgram : public Program {
    
    string output_filename;
    unsigned int iterations;

 public :
    TestProgram(char* ped, char* map, char* dat, char* outputfile, unsigned int iterations, bool verbose) : 
        Program(ped, map, dat, verbose), 
        output_filename(outputfile),
        iterations(iterations) {}
    
	~TestProgram() {}
    
    bool run();
};

#endif

