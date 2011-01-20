using namespace std;

#include <cstdio>

#include "lkg.h"
#include "person.h"
#include "pedigree.h"
#include "genotype.h"

bool Person::mendelian_errors() const {
	if(isfounder())
		return false;
		
	Person* m = ped->get_by_index(maternal_id);
	Person* p = ped->get_by_index(paternal_id);
		
	for(uint i = 0; i < genotypes.size(); ++i) {
		if(not genotype_compatible(m->get_genotype(i), 
								   p->get_genotype(i), 
									  get_genotype(i))) {
			fprintf(stderr, "error: %s, genotypes at loci number %d of \
				person \"%s\" inconsistent with parents\n", \
           		__func__, i+1, id.c_str());
			return true;
		}
	}
	
	return false;
}

void Person::fill_in_relationships() {
	Person* p;
	
	for(uint i = 0; i < ped->num_members(); ++i) {
		p = ped->get_by_index(i);
		
		if(p->get_mother() == id) {
			children.push_back(p);
			mates.push_back(ped->get_by_name(p->get_father()));
		}

		if(p->get_father() == id) {
			children.push_back(p);
			mates.push_back(ped->get_by_name(p->get_mother()));
		}
	}
}

string Person::gender_str() const {
	switch(gender) {
		case MALE :
			return "male";
		case FEMALE :
			return "female";
		case UNSEXED :
		default :
			break;
	}
	return "unspecified";
}

string Person::affection_str() const {
	switch(affection) {
		case UNAFFECTED :
			return "unaffected";
		case AFFECTED :
			return "affected";
		case UNKNOWN_AFFECTION :
		default :
			break;
	}
	return "unspecified";
}

void Person::print() const {
	printf("\tid: %s\n\
\tfather: %s\n\
\tmother: %s\n\
\tgender: %s\n\
\taffection: %s\n\
\tnumber of markers: %d\n",
			id.c_str(), 
			father.c_str(),
			mother.c_str(),
			gender_str().c_str(),
			affection_str().c_str(),
            int(genotypes.size()));
}

