#ifndef LKG_PEDIGREEPARSER_H_
#define LKG_PEDIGREEPARSER_H_

using namespace std;

#include <string>
#include <vector>

#include "genotype.h"
#include "person.h"
#include "parser.h"
#include "disease_model.h"


class Pedigree;
class GeneticMap;

class PedigreeParser : public Parser {

	vector<Pedigree>& pedigrees;
    DiseaseModel& disease_model;
	GeneticMap& map;

    bool ignore_genotypes;

	bool _parse_sex(const string& str, enum sex& s);
	bool _parse_affection(const string& str, enum affection& a);
	bool _parse_genotype(const string& a1, const string& a2, enum unphased_genotype& g);
	Pedigree& _current_ped(const string& famid);
	
 public :
	PedigreeParser(const string fn, vector<Pedigree>& peds, DiseaseModel& dm, GeneticMap& map) : 
	    Parser(fn, false), 
	    pedigrees(peds), 
	    disease_model(dm),
        map(map),
        ignore_genotypes(false) {}
    
    void set_ignore_genotypes() {
        ignore_genotypes = true;
    }
	bool parse_line(const int linenum, const string line);
    bool parse_end();
};

#endif

