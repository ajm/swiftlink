using namespace std;

#include <cstdio>
#include <vector>
#include <algorithm>
#include <queue>
#include <string>
#include <iostream>
#include <sstream>

#include "types.h"
#include "person.h"
#include "pedigree.h"
#include "peeling.h"
#include "disease_model.h"


Person::Person(const string name, const string father_name, const string mother_name, 
			   enum sex s, enum affection a, Pedigree* pedigree, const DiseaseModel& dm) :
        id(name),
    	mother(mother_name),
    	father(father_name),		
    	gender(s),
    	affection(a),			
    	internal_id(UNKNOWN_ID),
    	maternal_id(UNKNOWN_ID),
    	paternal_id(UNKNOWN_ID),
    	typed(false),
    	ped(pedigree),
    	genotypes(),
    	children(),
    	mates() {
    
    _init_probs(dm);
}
	
Person::Person(const Person& p) :
	    id(p.id),
		mother(p.mother),
		father(p.father),		
		gender(p.gender),
		affection(p.affection),			
		internal_id(p.internal_id),
		maternal_id(p.maternal_id),
		paternal_id(p.paternal_id),
		typed(p.typed),
		ped(p.ped),
		genotypes(p.genotypes),
		children(p.children),
		mates(p.mates) {
	    
	copy(p.disease_prob, p.disease_prob + 4, disease_prob);
}
	
Person& Person::operator=(const Person& rhs) {
	
	if(this != &rhs) {
	    id = rhs.id;
	    mother = rhs.mother;
	    father = rhs.father;
	    gender = rhs.gender;
	    affection = rhs.affection;
	    internal_id = rhs.internal_id;
	    maternal_id = rhs.maternal_id;
	    paternal_id = rhs.paternal_id;
	    typed = rhs.typed;
	    ped = rhs.ped;
	    genotypes = rhs.genotypes;
	    children = rhs.children;
	    mates = rhs.mates;
	    
	    copy(rhs.disease_prob, rhs.disease_prob + 4, disease_prob);
	}
	
	return *this;
}

void Person::_init_probs(const DiseaseModel& dm) {

    disease_prob[TRAIT_AA] = isfounder_str() ? \
        dm.get_apriori_prob(get_affection(), TRAIT_HOMO_A) : \
        dm.get_penetrance_prob(get_affection(), TRAIT_HOMO_A);

    disease_prob[TRAIT_AU] = \
    disease_prob[TRAIT_UA] = isfounder_str() ? \
        dm.get_apriori_prob(get_affection(), TRAIT_HETERO) : \
        dm.get_penetrance_prob(get_affection(), TRAIT_HETERO);

    disease_prob[TRAIT_UU] = isfounder_str() ? \
        dm.get_apriori_prob(get_affection(), TRAIT_HOMO_U) : \
        dm.get_penetrance_prob(get_affection(), TRAIT_HOMO_U);
}

bool Person::mendelian_errors() const {
	if(isfounder_str()) {
		return false;
    }
		
	Person* m = ped->get_by_index(maternal_id);
	Person* p = ped->get_by_index(paternal_id);
    
	for(unsigned int i = 0; i < genotypes.size(); ++i) {
        if(not genotype_compatible(m->get_genotype(i), 
								   p->get_genotype(i), 
									  get_genotype(i))) {
            
			fprintf(stderr, "error: genotypes at loci number %d of person \"%s\" inconsistent with parents\n", i+1, id.c_str());
			return true;
		}
	}
    
	return false;
}

void Person::add_mate(Person* p) {
    if(find(mates.begin(), mates.end(), p) == mates.end()) {
        mates.push_back(p);
    }
}

void Person::fill_in_relationships() {
	Person* p;
	
	for(unsigned int i = 0; i < ped->num_members(); ++i) {
		p = ped->get_by_index(i);
		
		if(p->get_mother() == id) {
			children.push_back(p);
			add_mate(ped->get_by_name(p->get_father()));
		}

		if(p->get_father() == id) {
			children.push_back(p);
		    add_mate(ped->get_by_name(p->get_mother()));
		}
	}
}

bool Person::is_offspring(unsigned int node) {
    for(unsigned i = 0; i < children.size(); ++i) {
        if(children[i]->get_internalid() == node) {
            return true;
        }
    }
    return false;
}

unsigned Person::count_unpeeled(vector<Person*>& v, PeelingState& ps) {
    unsigned int count = 0;

    for(unsigned int i = 0; i < v.size(); ++i) {
        if(not ps.is_peeled(v[i]->internal_id))
            ++count;
    }

    return count;
}

bool Person::offspring_peeled(PeelingState& ps) {
    return count_unpeeled(children, ps) == 0;
}

bool Person::partners_peeled(PeelingState& ps) {
    return count_unpeeled(mates, ps) == 0;
}

bool Person::safe_to_ignore_meiosis(enum parentage p) {
    Person* tmp = ped->get_by_index(p == MATERNAL ? maternal_id : paternal_id);
    
    if(not tmp->isfounder())
        return false;
        
    return tmp->num_children() == 1;
}

string Person::debug_string() {
    stringstream ss;
    
    //ss << "1\t" << internal_id + 1 << "\t" << paternal_id + 1 << "\t" << maternal_id + 1 << "\t" << gender << "\t" << affection;
    
    //return ss.str(); 
    
    ss.precision(DEBUG_FP_PRECISION);
    
    ss  << "\tid: " << id << "(" << internal_id << ")" << "\n" \
        << "\tfather: " << father << "(" << paternal_id << ")" << "\n" \
        << "\tmother: " << mother << "(" << maternal_id << ")" << "\n" \
        << "\tgender: " << gender_str(gender) << "\n" \
        << "\taffection: " << affection_str(affection) << "\n";
        
    //ss  << "\ttyped: " << typed << "\n";
    
    if(typed)
        ss << "\ttyped: yes\n";
    else
        ss << "\ttyped: no\n";
        //<< "\tnumber of markers: " << genotypes.size() << "\n";
    
    ss << "\tchildren:";
    for(unsigned i = 0; i < children.size(); ++i)
        ss << children[i]->get_internalid() << " ";
    ss << "\n";
    
    ss << "\tmates:";
    for(unsigned i = 0; i < mates.size(); ++i)
        ss << mates[i]->get_internalid() << " ";
    ss << "\n";
    
    /*
    ss << "\tprobabilities:" << "\n" \
       << "\t\tTRAIT_AA = " << fixed << disease_prob[TRAIT_AA] << "\n" \
       << "\t\tTRAIT_AU = " << fixed << disease_prob[TRAIT_AU] << "\n" \
       << "\t\tTRAIT_UA = " << fixed << disease_prob[TRAIT_UA] << "\n" \
       << "\t\tTRAIT_UU = " << fixed << disease_prob[TRAIT_UU] << "\n";
    */
    
    return ss.str(); 
}

