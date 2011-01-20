#ifndef LKG_LINKAGEPROGRAM_H_
#define LKG_LINKAGEPROGRAM_H_

#include "pedigree.h"
#include "geneticmap.h"
#include "diseasemodel.h"

class LinkageProgram {
	
	char *pedfile, 
		*mapfile, 
		*datfile;

	vector<Pedigree*> pedigrees;
	GeneticMap map;
	DiseaseModel dm;

	bool read_and_check_input();
	
 public :
	LinkageProgram(char* ped, char* map, char* dat)
		: pedfile(ped), mapfile(map), datfile(dat) {}
	~LinkageProgram();
    
	bool run();
	
	bool run_ga();
	bool run_sa();
};

#endif

