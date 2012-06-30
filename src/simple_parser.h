#ifndef LKG_SIMPLEPARSER_H_
#define LKG_SIMPLEPARSER_H_

#include <vector>
#include <string>
#include <cstdio>

#include "parser.h"


class SimpleParser : public Parser {

    vector<unsigned int> data;

 public :
	SimpleParser(const string fn) : 
	    Parser(fn, false),
		data() {}
    
	bool parse_line(const int linenum, const string s) {
		tokenise(s);
		
		if(tokens.size() == 0)
		    return true;
		
		if(tokens.size() != 1) {
		    fprintf(stderr, 
		            "error: read more than one field on line %d in file '%s' (read \"%s\")\n", 
		            linenum, filename.c_str(), s.c_str());
		    return false;
		}
		
		unsigned int tmp;
		stringstream ss(tokens[0]);
		
		if((ss >> tmp).fail()) {
		    fprintf(stderr, 
		            "error: not an integer, line %d, file '%s' (read \"%s\")\n", 
		            linenum, filename.c_str(), s.c_str());
		    return false;
		}
		
		data.push_back(tmp);
		
		return true;
	}
	
	vector<unsigned int>& get_values() {
	    return data;
	}
};

#endif

