using namespace std;

#include <cstdio>
#include <string>
#include <sstream>

#include "map_parser.h"
#include "genetic_map.h"


bool MapParser::parse_line(const int linenum, const string line) {
    string name;
    double gdist;
    int pdist;
    
    tokenise(line);
    
    if(tokens.size() == 0)
        return true;
    
    if(tokens.size() != 5) {
        fprintf(stderr, "error: %s, line %d: not enough data fields specified "
				"(expected at least 4, read %d)", 
				filename.c_str(), linenum+1, int(tokens.size()));
        return false;
    }
    
    // 0: chromosome, 1: genetic pos, 2: marker name, 3: physical pos, 4: index
    
    stringstream ss1(tokens[1]);
    if((ss1 >> gdist).fail()) {
        fprintf(stderr, "error: %s, line %d: genetic distance is not a floating"
				" point number (read '%s')\n", 
				filename.c_str(), linenum+1, tokens[1].c_str());
        return false;
    }
    
    if(gdist < 0.0) {
        fprintf(stderr, "error: %s, line %d: illegal genetic distance (%f)\n",
                filename.c_str(), linenum+1, gdist);
        return false;
    }
    
    name = tokens[2];
    
    stringstream ss2(tokens[3]);
    if((ss2 >> pdist).fail()) {
        fprintf(stderr, "error: %s, line %d: physical distance is not an "
				"integer (read '%s')\n", 
				filename.c_str(), linenum+1, tokens[3].c_str());
        return false;
    }
    
    if(pdist < 0) {
        fprintf(stderr, "error: %s, line %d: illegal physical distance (%d)\n",
                filename.c_str(), linenum+1, pdist);
        return false;
    }
    
    Snp snp(name, gdist, pdist);
    map.add(snp);
    
    return true;
}

