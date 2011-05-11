using namespace std;

#include <cstdio>
#include <string>
#include <fstream>

#include "linkage_writer.h"
#include "genetic_map.h"
#include "peeler.h"


bool LinkageWriter::write() {
	fstream hap;
	
	hap.open(filename.c_str(), ios::out);
    
	if(not hap.is_open()) {
		fprintf(stderr, "error: could not open linkage output file \"%s\"\n", 
			filename.c_str());
		return false;
	}
	
	for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
	    hap << i << "\t" << peeler->get(i) << endl;
	    
	    if(verbose) {
	        fprintf(stderr, "%d\t%.4f\n", i, peeler->get(i));
	    }
	}
	
	hap.close();
	
	return true;
}

