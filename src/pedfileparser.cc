using namespace std;

#include <cstdio>
#include <string>
#include <vector>

#include "pedfileparser.h"
#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "lkg.h"


bool PedigreeParser::_parse_sex(const string& str, enum sex* s) {
	if(str == "0")
		*s = UNSEXED;
	else if(str == "1")
		*s = MALE;
	else if(str == "2")
		*s = FEMALE;
	else
		return false;

	return true;
}

bool PedigreeParser::_parse_affection(const string& str, enum affection *a) {
	if(str == "0")
		*a = UNKNOWN_AFFECTION;
	else if(str == "1")
		*a = UNAFFECTED;
	else if(str == "2")
		*a = AFFECTED;
	else
		return false;

	return true;
}

bool PedigreeParser::_parse_genotype(const string& a1, const string& a2, 
		unphased_genotype_t *g) {
	if((a1 == "0") || (a2 == "0")) {
		*g = UNTYPED;
		return true;
	}
	
	if(a1 == a2) {
		if(a1 == "1") {
			*g = HOMOZ_A;
			return true;
		}
		if(a1 == "2") {
			*g = HOMOZ_B;
			return true;
		}
		
		return false;
	}
	else if(((a1 == "1") && (a2 == "2")) || \
			((a1 == "2") && (a2 == "1"))) {
		*g = HETERO;
		return true;
	}
	
	return false;
}

Pedigree* PedigreeParser::_current_ped(const string& famid) {
	Pedigree* p = NULL;

	for(unsigned int i = 0; i < pedigrees->size(); ++i) {
		if((*pedigrees)[i]->get_id() == famid) {
			//printf("\tfound family %s\n", famid.c_str());
			return (*pedigrees)[i];
		}
	}

	//printf("\tcreating new family %s\n", famid.c_str());
	p = new Pedigree(famid);
	pedigrees->push_back(p);

	return p;
}

bool PedigreeParser::parse_line(const int linenum, const string line) {
	string famid, id, fatherid, motherid;
	enum sex s = UNSEXED;
	enum affection a = UNKNOWN_AFFECTION;
	unphased_genotype_t g;
	Person* p = NULL;
	Pedigree* q = NULL;

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
				q = _current_ped(famid);
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
				if(!_parse_sex(tokens[i], &s)) {
					fprintf(stderr, "error: %s, line %d: bad sex \"%s\" (column 5)\n",
						filename.c_str(), linenum+1, tokens[i].c_str());
					return false;
				}
				break;
			case 5:
				if(!_parse_affection(tokens[i], &a)) {
					fprintf(stderr, "error: %s, line %d: bad affection status \
\"%s\" (column 6)\n", filename.c_str(), linenum+1, tokens[i].c_str());
					return false;
				}
				break;
			default:
				break;
		}
	}
	
	p = new Person(id, fatherid, motherid, s, a, q);
//	if(! q->add(p)) {
//		delete p;
//		return false;
//	}
	
	// handle genotypes
	for(int i = 0; i < ((tokens.size() - 6) / 2.0); ++i) {
		int index = 6 + (i * 2);
		if(!_parse_genotype(tokens[index], tokens[index+1], &g)) {
			fprintf(stderr, "error: %s, line %d: error parsing genotype %d \
(columns %d and %d)\n", filename.c_str(), linenum+1, i+1, index+1, index+2);
			return false;
		}
		p->add_genotype(g);
	}
    
    if(! q->add(p)) {
		delete p;
		return false;
	}    

    delete p;
	return true;
}

