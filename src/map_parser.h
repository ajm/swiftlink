#ifndef LKG_MAPPARSER_H_
#define LKG_MAPPARSER_H_

using namespace std;

#include <cstdio>
#include <string>
#include <vector>

#include "parser.h"
#include "genetic_map.h"


class MapParser : public Parser {
    
    GeneticMap& map;
    
 public :
    MapParser(const string fn, GeneticMap& m) 
        : Parser(fn, false), map(m) {}
    
    bool parse_line(const int linenum, const string line);
    bool parse_end();
};

#endif

