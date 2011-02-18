#ifndef LKG_PROGRAM_H_
#define LKG_PROGRAM_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "pedigree.h"
#include "geneticmap.h"
#include "disease_model.h"


class Program {

 private :
    char* pedfile;
    char* mapfile;
	char* datfile;

	bool read_and_check_input();

 protected :
	vector<Pedigree> pedigrees;
	GeneticMap map;
	DiseaseModel dm;
	
 public :
	Program(char* ped, char* map, char* dat)
		: pedfile(ped), mapfile(map), datfile(dat) {

        if(! read_and_check_input()) {
            exit(1);
        }
    }
    
	virtual ~Program() {}
    
	virtual bool run()=0;
};

#endif

