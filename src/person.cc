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
#include "genetic_map.h"

using namespace std;


Person::Person(const string name, const string father_name, const string mother_name, 
			   enum sex s, enum affection a, Pedigree* pedigree, DiseaseModel& dm) :
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
        dm(&dm),
        genotypes(),
        genotypes_prob(),
        children(),
        mates() {
    
    _init_probs();
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
        dm(p.dm),
		genotypes(p.genotypes),
        genotypes_prob(p.genotypes_prob),
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
        dm = rhs.dm;
	    genotypes = rhs.genotypes;
        genotypes_prob = rhs.genotypes_prob;
	    children = rhs.children;
	    mates = rhs.mates;
	    
	    copy(rhs.disease_prob, rhs.disease_prob + 4, disease_prob);
	}
	
	return *this;
}

void Person::_init_probs() {
/*
    disease_prob[TRAIT_AA] = isfounder_str() ? \
        dm->get_apriori_prob(get_affection(), TRAIT_HOMO_A) : \
        dm->get_penetrance_prob(get_affection(), TRAIT_HOMO_A);

    disease_prob[TRAIT_AU] = \
    disease_prob[TRAIT_UA] = isfounder_str() ? \
        dm->get_apriori_prob(get_affection(), TRAIT_HETERO) : \
        dm->get_penetrance_prob(get_affection(), TRAIT_HETERO);

    disease_prob[TRAIT_UU] = isfounder_str() ? \
        dm->get_apriori_prob(get_affection(), TRAIT_HOMO_U) : \
        dm->get_penetrance_prob(get_affection(), TRAIT_HOMO_U);

    if(dm->is_sexlinked() and ismale()) {
        disease_prob[TRAIT_AU] = 0.0;
        disease_prob[TRAIT_UA] = 0.0;
    }
*/
    disease_prob[TRAIT_AA] = isfounder_str() ? \
        dm->get_apriori_prob2(get_affection(), TRAIT_HOMO_A, get_sex()) : \
        dm->get_penetrance_prob2(get_affection(), TRAIT_HOMO_A, get_sex());

    disease_prob[TRAIT_AU] = \
    disease_prob[TRAIT_UA] = isfounder_str() ? \
        dm->get_apriori_prob2(get_affection(), TRAIT_HETERO, get_sex()) : \
        dm->get_penetrance_prob2(get_affection(), TRAIT_HETERO, get_sex());

    disease_prob[TRAIT_UU] = isfounder_str() ? \
        dm->get_apriori_prob2(get_affection(), TRAIT_HOMO_U, get_sex()) : \
        dm->get_penetrance_prob2(get_affection(), TRAIT_HOMO_U, get_sex());

    //printf("XXX %s %f %f %f %f\n", id.c_str(), disease_prob[TRAIT_AA], disease_prob[TRAIT_AU], disease_prob[TRAIT_UA], disease_prob[TRAIT_UU]);
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
									  get_genotype(i),
                                      get_sex(),
                                  dm->is_sexlinked())) {
            
			fprintf(stderr, "error: genotypes at loci number %d of person \"%s\" inconsistent with parents\n", i+1, id.c_str());
			return true;
		}
	}
    
	return false;
}

/*
bool Person::all_homozygous() {

    for(unsigned int i = 0; i < genotypes.size(); ++i) {
        if(get_genotype(i) == HETERO)
            return false;
    }

    return true;
}
*/

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
    
    if(not tmp->isfounder()) {
        if(dm->is_sexlinked()) {
            // sex == MALE   : paternal Y, maternal can be either
            // sex == FEMALE : paternal grandmothers allele, maternal can be either
            return p == PATERNAL;
        }

        return false;
    }
        
    return tmp->num_children() == 1;
}

void Person::populate_trait_prob_cache(GeneticMap& map) {
    vector<double> probs(4, 0.0);

    genotypes_prob.clear();
    
    for(unsigned int i = 0; i < genotypes.size(); ++i) {
        for(int j = 0; j < 4; ++j) {
            enum phased_trait pt = static_cast<enum phased_trait>(j);
            double marker_prob = map.get_prob(i, pt, dm->is_sexlinked() and ismale());

            probs[j] = get_genotype_probability(genotypes[i], pt, marker_prob);
        }

        double total = probs[0] + probs[1] + probs[2] + probs[3];

        for(int j = 0; j < 4; ++j) {
            probs[j] /= total;
        }

        genotypes_prob.push_back(probs);
    }
}

double Person::get_genotype_probability(enum unphased_genotype g, enum phased_trait pt, double marker_prob) {
    
    // penetrace prob
    if(not isfounder()) {
        if(istyped()) {
            switch(g) {
                case HETERO :
                    return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? 1.0 : 0.0;
                case HOMOZ_A :
                    return (pt == TRAIT_UU) ? 1.0 : 0.0;
                case HOMOZ_B :
                    return (pt == TRAIT_AA) ? 1.0 : 0.0;
                default :
                    return 1.0;
            }
        }
        else {
            if(dm->is_sexlinked()) {
                if(ismale()) {
                    return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? 0.0 : 1.0;
                }
            }

            return 1.0;
        }
    }
    
    // penetrance + founder prob
    if(istyped()) {
        switch(g) {
            case HETERO :
                return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? marker_prob : 0.0;
            case HOMOZ_A :
                return (pt == TRAIT_UU) ? marker_prob : 0.0;
            case HOMOZ_B :
                return (pt == TRAIT_AA) ? marker_prob : 0.0;
            default :
                return marker_prob;
        }
    }
    else {
        if(dm->is_sexlinked()) {
            if(ismale()) {
                return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? 0.0 : marker_prob;
            }
        }

        return marker_prob;
    }
    
    fprintf(stderr, "error: %s:%d\n", __FILE__, __LINE__);
    abort();
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

