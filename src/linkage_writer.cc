using namespace std;

#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>

#include "linkage_writer.h"
#include "genetic_map.h"
#include "peeler.h"


bool LinkageWriter::write(vector<Peeler*>& peelers) {
    fstream f;
    	
	f.open(filename.c_str(), fstream::out | fstream::trunc);
    
	if(not f.is_open()) {
		fprintf(stderr, "error: could not open linkage output file \"%s\"\n", filename.c_str());
		return false;
	}
	
	for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
        double tmp = peelers[0]->get(i);
        for(unsigned j = 1; j < peelers.size(); ++j) {
            tmp += peelers[j]->get(i);
        }
        
        stringstream ss;
        
        // marker id + total LOD
	    ss << map->get_name(i) << "\t" << tmp;
        
        // more than one family
        if(peelers.size() > 1) {
            for(unsigned j = 0; j < peelers.size(); ++j) {
                ss << "\t" << peelers[j]->get(i);
            }
        }
        
        ss << "\n";
        f << ss.str();
        
        if(verbose)
            fprintf(stderr, "%s", ss.str().c_str());
	}
	
	f.close();
	
	return true;
}
