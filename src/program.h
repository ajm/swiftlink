#ifndef LKG_PROGRAM_H_
#define LKG_PROGRAM_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "pedigree.h"
#include "genetic_map.h"
#include "disease_model.h"
#include "types.h"


class Program {

 private :
    string pedfile;
    string mapfile;
	string datfile;

	bool read_and_check_input();

 protected :
	vector<Pedigree> pedigrees;
	GeneticMap map;
	DiseaseModel dm;
	struct mcmc_options options;
	string outfile;
	
 public :
	Program(const char* ped, const char* map, const char* dat, const char* out, struct mcmc_options options) : 
	    pedfile(ped), 
	    mapfile(map), 
	    datfile(dat),
	    pedigrees(),
	    map(options.lodscores),
	    dm(),
	    options(options),
	    outfile(out) {

        if(! read_and_check_input()) {
            exit(1);
        }
    }
    
    Program(const Program& rhs) :
        pedfile(rhs.pedfile),
        mapfile(rhs.mapfile), 
	    datfile(rhs.datfile), 
	    pedigrees(rhs.pedigrees),
	    map(rhs.map),
	    dm(rhs.dm),
	    options(rhs.options),
	    outfile(rhs.outfile) {}
    
    Program& operator=(const Program& rhs) {
        
        if(&rhs != this) {
            pedfile = rhs.pedfile;
            mapfile = rhs.mapfile; 
	        datfile = rhs.datfile; 
	        pedigrees = rhs.pedigrees;
	        map = rhs.map;
	        dm = rhs.dm;
	        options = rhs.options;
	        outfile = rhs.outfile;
        }
        
        return *this;
    }
    
	virtual ~Program() {}
	
	virtual bool run()=0;
};

#endif

