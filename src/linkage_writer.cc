using namespace std;

#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>

#include "linkage_writer.h"
#include "genetic_map.h"
#include "peeler.h"


bool LinkageWriter::write(vector<double*>& lod_scores) {
    fstream f;
    	
	f.open(filename.c_str(), fstream::out | fstream::trunc);
    
	if(not f.is_open()) {
		fprintf(stderr, "error: could not open linkage output file \"%s\"\n", filename.c_str());
		return false;
	}
	
	for(unsigned int i = 0; i < (map->num_markers() - 1); ++i) {
        double tmp = lod_scores[0][i];
        for(unsigned int j = 1; j < lod_scores.size(); ++j) {
            tmp += lod_scores[j][i];
        }
        
        stringstream ss;
        
        // marker id + total LOD
	    ss << map->get_name(i) << "\t" << tmp;
        
        // more than one family
        if(lod_scores.size() > 1) {
            for(unsigned int j = 0; j < lod_scores.size(); ++j) {
                ss << "\t" << lod_scores[j][i];
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
