#ifndef LKG_LINKAGEPARSER_H_
#define LKG_LINKAGEPARSER_H_

using namespace std;

#include <string>
#include <vector>

#include "parser.h"

class GeneticMap;
class DiseaseModel;

class LinkageParser : public Parser {
	
    GeneticMap* map;
	DiseaseModel* dis;

	int marker_linenum;
	int recomb_linenum;
	int marker_code;
	double marker_freq;
	int marker_alleles;
	int markers_read[4];
	int liability_classes;
	
	int number_of_loci;
	int program_code;
		
	bool general_information(const int linenum);
	bool description_of_loci(const int linenum);
	bool information_on_recombination(const int linenum);
	bool set_numloci(const string s);
	bool set_sexlinked(const string s);
	bool set_programcode(const string s);
	bool ensure_false(const string s, const string msg);
	bool check_num_loci();
	bool read_affection_status();
	bool read_numbered_allele();
	bool read_binary_factor();
	bool read_quantitative_variable();
	bool unsupported(const string s);
	bool read_marker_code(const string s);
	bool read_abstract_marker();
	vector<double>* get_doubles();
	int total_markers_read();
	void marker_end();

 public :
    LinkageParser(const string fn, GeneticMap* m, DiseaseModel* d) 
        : Parser(fn, true), map(m), dis(d), marker_linenum(0), recomb_linenum(0), 
		  marker_code(-1), number_of_loci(-1), program_code(-1) {
		for(int i = 0; i < 4; ++i)
			markers_read[i] = 0;
	}
	
    bool parse_line(const int linenum, const string line);
	bool parse_end();
};


#endif

