#ifndef LKG_CONFIGPARSER_H_
#define LKG_CONFIGPARSER_H_

using namespace std;

#include <cstdio>
#include <string>
#include <sstream>

class ConfigParser : public Parser {

 public :
    unsigned population;
    unsigned iterations;
    unsigned elitism;
    
    int crossover;

    int selection;
    unsigned selection_size;
    double selection_prob;
    
 
	ConfigParser(const string fn) 
		: Parser(fn, true) {}

	bool parse_line(const int linenum, const string s) {
	    stringstream ss;
	
		tokenise(s);
		
		if(tokens[0] == "population") {
		    ss << tokens[1];
		    if((ss >> population).fail()) {
		        fprintf(stderr, "bad config %d (%s)\n", linenum, s.c_str());
		        return false;
		    }
		}
		
		if(tokens[0] == "iterations") {
		    ss << tokens[1];
		    if((ss >> iterations).fail()) {
		        fprintf(stderr, "bad config %d (%s)\n", linenum, s.c_str());
		        return false;
		    }
		}
		
		if(tokens[0] == "elitism") {
		    ss << tokens[1];
		    if((ss >> elitism).fail()) {
		        fprintf(stderr, "bad config %d (%s)\n", linenum, s.c_str());
		        return false;
		    }
		}
		
		if(tokens[0] == "crossover") {
		    if(tokens[1] == "single") {
		        crossover = 0;
		        return true;
		    }
		    if(tokens[1] == "double") {
		        crossover = 1;
		        return true;
		    }
		    if(tokens[1] == "uniform") {
		        crossover = 2;
		        return true;
		    }
		    fprintf(stderr, "bad config %d (%s)\n", linenum, s.c_str());
		    return false;
		}
		
		if(tokens[0] == "selection") {
		    if(tokens[1] == "truncation") {
		        selection = 0;
		        ss << tokens[2];
		        if((ss >> selection_prob).fail()) {
		            fprintf(stderr, "bad config %d (%s)\n", linenum, s.c_str());
		            return false;
		        }
		    }
		    if(tokens[1] == "tournament") {
		        selection = 1;
		        ss << tokens[2];
		        if((ss >> selection_size).fail()) {
		            fprintf(stderr, "bad config %d (%s)\n", linenum, s.c_str());
		            return false;
		        }
		        ss.clear();
		        ss << tokens[3];
		        if((ss >> selection_prob).fail()) {
		            fprintf(stderr, "bad config %d (%s) %s\n", linenum, s.c_str(), tokens[3].c_str());
		            return false;
		        }
		    }
		}
		
		return true;
	}
};

#endif

