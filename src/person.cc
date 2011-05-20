using namespace std;

#include <cstdio>
#include <vector>
#include <algorithm>
#include <queue>
#include <string>
#include <iostream>
#include <sstream>

#include "lkg.h"
#include "person.h"
#include "pedigree.h"
#include "genotype.h"
#include "peeling.h"
#include "disease_model.h"
#include "trait.h"
#include "debug.h"


void Person::init_probs(DiseaseModel& dm) {

    disease_prob[TRAIT_AA] = isfounder() ? \
        dm.get_apriori_prob(get_affection(), TRAIT_HOMO_A) : \
        dm.get_penetrance_prob(get_affection(), TRAIT_HOMO_A);

    disease_prob[TRAIT_AU] = \
    disease_prob[TRAIT_UA] = isfounder() ? \
        dm.get_apriori_prob(get_affection(), TRAIT_HETERO) : \
        dm.get_penetrance_prob(get_affection(), TRAIT_HETERO);

    disease_prob[TRAIT_UU] = isfounder() ? \
        dm.get_apriori_prob(get_affection(), TRAIT_HOMO_U) : \
        dm.get_penetrance_prob(get_affection(), TRAIT_HOMO_U);
}

double Person::get_disease_prob(enum phased_trait pt) {
    return disease_prob[pt];
}

bool Person::mendelian_errors() const {
	if(isfounder()) {
		return false;
    }
		
	Person* m = ped->get_by_index(maternal_id);
	Person* p = ped->get_by_index(paternal_id);
		
	for(unsigned int i = 0; i < genotypes.size(); ++i) {
		if(not genotype_compatible(m->get_genotype(i), 
								   p->get_genotype(i), 
									  get_genotype(i))) {
			fprintf(stderr, "error: %s, genotypes at loci number %d of person \"%s\" inconsistent with parents\n",
           		__func__, i+1, id.c_str());
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
			//mates.push_back(ped->get_by_name(p->get_father()));
			add_mate(ped->get_by_name(p->get_father()));
		}

		if(p->get_father() == id) {
			children.push_back(p);
			//mates.push_back(ped->get_by_name(p->get_mother()));
		    add_mate(ped->get_by_name(p->get_mother()));
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

string Person::debug_string() {
    stringstream ss;
    
    ss.precision(DEBUG_FP_PRECISION);
    
    ss << "\tid: " << id << "(" << internal_id << ")" << endl \
       << "\tfather: " << father << "(" << paternal_id << ")" << endl \
       << "\tmother: " << mother << "(" << maternal_id << ")" << endl \
       << "\tgender: " << gender_str() << endl \
       << "\taffection: " << affection_str() << endl \
       << "\tnumber of markers: " << genotypes.size() << endl;
       
    ss << "\tchildren:";
    for(unsigned i = 0; i < children.size(); ++i)
        ss << children[i]->internal_id << " ";
    ss << endl;
    
    ss << "\tmates:";
    for(unsigned i = 0; i < mates.size(); ++i)
        ss << mates[i]->internal_id << " ";
    ss << endl;
    
    ss << "\tprobabilities:" << endl \
       << "\t\tTRAIT_AA = " << fixed << disease_prob[TRAIT_AA] << endl \
       << "\t\tTRAIT_AU = " << fixed << disease_prob[TRAIT_AU] << endl \
       << "\t\tTRAIT_UA = " << fixed << disease_prob[TRAIT_UA] << endl \
       << "\t\tTRAIT_UU = " << fixed << disease_prob[TRAIT_UU] << endl;
       
    return ss.str(); 
}

//----------
bool Person::contains_unpeeled(vector<Person*>& v, PeelingState& ps) {

    for(unsigned int i = 0; i < v.size(); ++i) {
        if(not ps.is_peeled(v[i]->internal_id))
            return true;
    }

    return false;
}

unsigned int Person::count_unpeeled(vector<Person*>& v, PeelingState& ps) {
    unsigned int count = 0;

    for(unsigned int i = 0; i < v.size(); ++i) {
        if(not ps.is_peeled(v[i]->internal_id))
            ++count;
    }

    return count;
}

bool Person::offspring_peeled(PeelingState& ps) {
    return not contains_unpeeled(children, ps);
}

bool Person::partners_peeled(PeelingState& ps) {
    return not contains_unpeeled(mates, ps);
}

bool Person::parents_peeled(PeelingState& ps) {
    return isfounder() or (ps.is_peeled(maternal_id) and ps.is_peeled(paternal_id));
}

bool Person::ripe_above(PeelingState& ps) {
    return parents_peeled(ps);
}

bool Person::ripe_below(PeelingState& ps) {
    return not isfounder() and offspring_peeled(ps);
}

bool Person::ripe_across(PeelingState& ps) {
    return  offspring_peeled(ps) and \
            parents_peeled(ps) and \
            (count_unpeeled(mates, ps) == 1);
}

bool Person::ripe_all(PeelingState& ps) {
    return  offspring_peeled(ps) and \
            parents_peeled(ps) and \
            partners_peeled(ps);
}

bool Person::ripe_above_partners(PeelingState& ps) {
    
    for(unsigned int i = 0; i < mates.size(); ++i) {
        if(not mates[i]->ripe_above(ps))
            return false;
    }

    return ripe_above(ps);
}

bool Person::ripe_to_peel_down(PeelingState& ps) {
    return (not ps.is_peeled(internal_id)) and \
            ripe_above(ps) and \
           (count_unpeeled(mates, ps) == 1) and \
           (count_unpeeled(children, ps) == 1);
}

bool Person::ripe_parents(PeelingState& ps) {

    if(isfounder())
        return false;
    
    Person* mother = ped->get_by_index(maternal_id);
    Person* father = ped->get_by_index(paternal_id);

    return mother->ripe_to_peel_down(ps) and father->ripe_to_peel_down(ps);   
}

bool Person::ripe(PeelingState& ps) {
    return ripe_all(ps) or \
            ripe_across(ps) or \
            ripe_below(ps) or \
            ripe_above_partners(ps);
}

unsigned int Person::get_unpeeled_spouse(PeelingState& ps) {
    vector<unsigned int> unpeeled_spouses;
    
    for(unsigned i = 0; i < mates.size(); ++i) {
        unsigned int spouse_id = mates[i]->get_internalid();
        if(not ps.is_peeled(spouse_id)) {
            unpeeled_spouses.push_back(spouse_id);
        }
    }
    
    if(unpeeled_spouses.size() != 1) {
        fprintf(stderr, "error: a PARENT_PEEL was considered legal, but there are %d unpeeled spouses\n", int(unpeeled_spouses.size()));
        abort();
    }
    
    return unpeeled_spouses[0];
}

bool Person::peel_operation(PeelOperation& p, PeelingState& state) {
    p.set_pivot(internal_id);
    p.set_type(NULL_PEEL);

/*
    // would be nice, but then I don't know what to do later...
    if(ripe(state)) {
        get_cutset(p, state);
        return true;
    }
*/
    
    if(ripe_all(state)) {
        p.set_type(LAST_PEEL);
    }
    else if(ripe_across(state)) {
        p.set_type(PARTNER_PEEL);
    }
    else if(ripe_below(state)) {
        p.set_type(CHILD_PEEL);
    }
/*
    else if(ripe_parents(state)) {
        p.set_type(PARENT_PEEL);
    }
*/
    if(p.get_type() != NULL_PEEL) {
        get_cutset(p, state);
        return true;
    }
    
    return false;
}

void Person::neighbours(vector<unsigned int>& nodes) {
    nodes.clear();

    if(maternal_id != UNKNOWN_ID)
        nodes.push_back(maternal_id);

    if(paternal_id != UNKNOWN_ID)
        nodes.push_back(paternal_id);

    for(unsigned int i = 0; i < children.size(); ++i)
        nodes.push_back(children[i]->internal_id);

    for(unsigned int i = 0; i < mates.size(); ++i)
        nodes.push_back(mates[i]->internal_id);
}

void Person::get_cutset(PeelOperation& operation, PeelingState& state) {
    queue<unsigned int> q;
    vector<unsigned int> n;
    //unsigned int visited[ped->num_members()]; // ISO C++ bitches...
    vector<int> visited(ped->num_members(), WHITE);
    unsigned int tmp, tmp2;
    Person* p;
/*
    for(unsigned int i = 0; i < ped->num_members(); ++i) {
		visited[i] = WHITE;
    }
*/    

    if(operation.get_type() == PARENT_PEEL) {
        visited[maternal_id] = GREY;
        visited[paternal_id] = GREY;
        
        q.push(maternal_id);
        q.push(paternal_id);
    }
    else {
        visited[internal_id] = GREY;
        q.push(internal_id);
    }

    while(not q.empty()) {
		tmp = q.front();
		q.pop();
		
        p = ped->get_by_index(tmp);
        p->neighbours(n);

        for(unsigned int i = 0; i < n.size(); ++i) {
            tmp2 = n[i];

            if(not state.is_peeled(tmp2)) {
                if((internal_id != tmp2) or (operation.get_type() == PARENT_PEEL)) {
                    operation.add_cutnode(tmp2);
                }
                continue;
            }

            if(visited[tmp2] == WHITE) {
                visited[tmp2] = GREY;
                q.push(tmp2);
            }
        }
        
        visited[tmp] = BLACK;
    }
    
    /*
    if(operation.get_type() == PARENT_PEEL) {
        operation.add_cutnode(maternal_id);
        operation.add_cutnode(paternal_id);
    }
    */
}

