using namespace std;

#include <cstdio>
#include <string>
#include <vector>

#include "pedigree_parser.h"
#include "pedigree.h"
#include "person.h"
#include "types.h"


bool PedigreeParser::_parse_sex(const string& str, enum sex& s) {
	if(str == "0")
		s = UNSEXED;
	else if(str == "1")
		s = MALE;
	else if(str == "2")
		s = FEMALE;
	else
		return false;

	return true;
}

bool PedigreeParser::_parse_affection(const string& str, enum affection &a) {
	if(str == "0")
		a = UNKNOWN_AFFECTION;
	else if(str == "1")
		a = UNAFFECTED;
	else if(str == "2")
		a = AFFECTED;
	else
		return false;

	return true;
}

bool PedigreeParser::_parse_genotype(const string& a1, const string& a2, enum unphased_genotype& g) {
	if((a1 == "0") || (a2 == "0")) {
		g = UNTYPED;
		return true;
	}
	
	if(a1 == a2) {
		if(a1 == "1") {
			g = HOMOZ_A;
			return true;
		}
		if(a1 == "2") {
			g = HOMOZ_B;
			return true;
		}
		
		return false;
	}
	else if(((a1 == "1") && (a2 == "2")) || \
			((a1 == "2") && (a2 == "1"))) {
		g = HETERO;
		return true;
	}
	
	return false;
}

Pedigree& PedigreeParser::_current_ped(const string& famid) {

	for(unsigned int i = 0; i < pedigrees.size(); ++i) {
		if(pedigrees[i].get_id() == famid) {
			return pedigrees[i];
		}
	}

    Pedigree p(famid);
	pedigrees.push_back(p);
        
	return pedigrees.back();
}

bool PedigreeParser::parse_line(const int linenum, const string line) {
	string famid, id, fatherid, motherid;
	enum sex s = UNSEXED;
	enum affection a = UNKNOWN_AFFECTION;
	enum unphased_genotype g;
	
	tokenise(line);

	//printf("read: %s\n", line.c_str());

	if(tokens.size() == 0) // empty line
		return true;

	if(tokens.size() < 6) {
		fprintf(stderr,"error: %s, line %d: not enough fields, only %d found\n",
					filename.c_str(), linenum+1, int(tokens.size()));
		return false;
	}

	if((tokens.size() % 2) != 0) {
		fprintf(stderr,"error: %s, line %d: contains an odd number of alleles\n",
					filename.c_str(), linenum+1);
		return false;
	}
	
	// info that defines person
	for(int i = 0; i < 6; ++i) {
		switch(i) {
			case 0:
				famid = tokens[i];
				break;
			case 1:
				id = tokens[i];
				break;
			case 2:
				fatherid = tokens[i];
				break;
			case 3:
				motherid = tokens[i];
				break;
			case 4:
				if(not _parse_sex(tokens[i], s)) {
					fprintf(stderr, "error: %s, line %d: bad sex \"%s\" (column 5)\n",
						filename.c_str(), linenum+1, tokens[i].c_str());
					return false;
				}
				break;
			case 5:
				if(not _parse_affection(tokens[i], a)) {
					fprintf(stderr, "error: %s, line %d: bad affection status \"%s\" (column 6)\n", 
					            filename.c_str(), linenum+1, tokens[i].c_str());
					return false;
				}
				break;
			default:
				break;
		}
	}

    // create / retrieve objects	
    Pedigree& q = _current_ped(famid);
	Person p(id, fatherid, motherid, s, a, &q, disease_model);
    	
    if(not ignore_genotypes) {
	    // handle genotypes
	    for(int i = 0; i < ((tokens.size() - 6) / 2.0); ++i) {
		    int index = 6 + (i * 2);

		    if(not _parse_genotype(tokens[index], tokens[index+1], g)) {
			    fprintf(stderr, "error: %s, line %d: error parsing genotype %d (columns %d and %d)\n", 
                    filename.c_str(), linenum+1, i+1, index+1, index+2);
			    return false;
		    }

		    p.add_genotype(g);
	    }

        p.populate_trait_prob_cache(map);
    }

    return q.add(p);
}

bool PedigreeParser::parse_end() {
    bool ret = true;

    if(pedigrees.size() == 0) {
        fprintf(stderr, "error: %s contains no families\n", 
            filename.c_str());
            
        return false;
    }
    
    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        
		if(not pedigrees[i].sanity_check()) {
            fprintf(stderr, "error: %s, pedigree %s contains relationship errors\n",
                filename.c_str(), pedigrees[i].get_id().c_str());
			ret = false;
        }
    }

    return ret;
}

