using namespace std;

#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>

#include "linkage_writer.h"
#include "genetic_map.h"
#include "peeler.h"
#include "lod_score.h"


bool LinkageWriter::write(vector<LODscores*>& all_scores) {
    fstream f;
    	
	f.open(filename.c_str(), fstream::out | fstream::trunc);
    
	if(not f.is_open()) {
		fprintf(stderr, "error: could not open linkage output file \"%s\"\n", filename.c_str());
		return false;
	}
	
	/*
	// old version
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
	*/
	
	for(unsigned int i = 0; i < (map->num_markers() - 1); ++i) {
	    for(unsigned int j = 0; j < map->get_lodscore_count(); ++j) {
	        double tmp = all_scores[0]->get(i, j);
            for(unsigned int k = 1; k < all_scores.size(); ++k) {
                tmp += all_scores[k]->get(i,j);
            }
            
            stringstream ss;
            
            // genetic position + total LOD
            // (get_genetic_position needs j+1 because it treats j as an offset, ie: indexed from 1)
            ss << 100 * map->get_genetic_position(i, j+1) << "\t" << tmp;
            
            // more than one family
            if(all_scores.size() > 1) {
                for(unsigned int k = 0; k < all_scores.size(); ++k) {
                    ss << "\t" << all_scores[k]->get(i,j);
                }
            }
            
            ss << "\n";
            f << ss.str();
            
            if(verbose)
                fprintf(stderr, "%s", ss.str().c_str());
	    }
	}
	
	f.close();
	
	return true;
}
