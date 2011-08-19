using namespace std;

#include <cstdio>
#include <vector>
#include <algorithm>
#include <queue>
#include <string>
#include <iostream>
#include <sstream>

#include "misc.h"
#include "person.h"
#include "pedigree.h"
#include "genotype.h"
#include "peeling.h"
#include "disease_model.h"
#include "trait.h"


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
    	typed(true),
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
	    
	copy(p.disease_prob, p.disease_prob+4, disease_prob);
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

double Person::get_disease_prob(enum phased_trait pt) {
    return disease_prob[pt];
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

bool Person::founder_mates_peeled(PeelingState& ps) {
    for(unsigned int i = 0; i < mates.size(); ++i) {
        if(mates[i]->isfounder() and (not ps.is_peeled(mates[i]->internal_id)) and (not is_parent(i)))
            return false;
    }

    return true;
}

bool Person::partners_peeled(PeelingState& ps) {
    return count_unpeeled(mates, ps) == 0;
}

bool Person::parents_peeled(PeelingState& ps) {
    return isfounder() or (ps.is_peeled(maternal_id) and ps.is_peeled(paternal_id));
}

bool Person::one_parent_peeled(PeelingState& ps) {
    return ps.is_peeled(maternal_id) xor ps.is_peeled(paternal_id);
}

bool Person::ripe_above(PeelingState& ps) {
    return parents_peeled(ps);
}

bool Person::ripe_above_at_least_one_parent(PeelingState& ps) {
    return ped->get_by_index(maternal_id)->ripe_above(ps) or \
           ped->get_by_index(paternal_id)->ripe_above(ps);
}

bool Person::ripe_to_peel_up(PeelingState& ps) {
    return not isfounder() and \
           not parents_peeled(ps) and \
           offspring_peeled(ps) and \
           founder_mates_peeled(ps);
}

bool Person::ripe_to_peel_across(PeelingState& ps) {
    return  (parents_peeled(ps) and offspring_peeled(ps) and (count_unpeeled(mates, ps) == 1)) or \
            (parents_peeled(ps) and partners_peeled(ps)) or \
            (not isfounder() and one_parent_peeled(ps) and offspring_peeled(ps) and founder_mates_peeled(ps));
}


bool Person::ripe_to_peel_final(PeelingState& ps) {
    //return offspring_peeled(ps) and parents_peeled(ps) and partners_peeled(ps);
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        if((i != internal_id) and not ps.is_peeled(i))
            return false;
    }
    
    return true;
}

bool Person::ripe_above_singular_mate(PeelingState& ps) {
    Person* p = NULL;
    int count = 0;
    
    for(unsigned i = 0; i < mates.size(); ++i) {
        if(not ps.is_peeled(mates[i]->get_internalid())) {
            p = mates[i];
            count++;
        }
    }
    
    if(count != 1)
        return false;
    
    return p->ripe_above(ps);
}

unsigned Person::get_unpeeled_mate(PeelingState& ps) {
    
    for(unsigned i = 0; i < mates.size(); ++i) {
        unsigned pid = mates[i]->get_internalid();
        
        if(not ps.is_peeled(pid)) {
            return pid;
        }
    }
    
    return UNKNOWN_ID;
}

/*
 
 it must be:
    ripe above this node
    there only be one unpeeled mate node
    ripe above the unpeeled mate node
    one child to peel down to (XXX temporary)
    that child is via the unpeeled mate
 
 */

bool Person::ripe_to_peel_down(PeelingState& ps) {
    return ripe_above(ps) and \
        ripe_above_singular_mate(ps);// and (count_unpeeled(children, ps) == 1);
}

bool Person::peel_operation(PeelOperation& p, PeelingState& state) {
    if(state.is_peeled(internal_id)) {
        return false;
    }
    
    //p.set_pivot(internal_id);
    p.set_type(NULL_PEEL);

/*
    // would be nice, but then I don't know what to do later...
    if(ripe(state)) {
        get_cutset(p, state);
        return true;
    }
*/
    
    if(ripe_to_peel_final(state)) {
        p.set_type(LAST_PEEL);
    }
    else if(ripe_to_peel_across(state)) {
        p.set_type(PARTNER_PEEL);
    }
    else if(ripe_to_peel_up(state)) {
        p.set_type(CHILD_PEEL);
    }
    
    else if(ripe_to_peel_down(state)) {
        p.set_type(PARENT_PEEL);
        //p.add_peelnode(get_unpeeled_mate(state));
    }
    
    if(p.get_type() != NULL_PEEL) {
        p.set_peelnode(internal_id);
        get_cutset(p, state);
        return true;
    }
    
    return false;
}

void Person::neighbours(vector<unsigned int>& nodes, PeelingState& ps) {
    nodes.clear();

    if(maternal_id != UNKNOWN_ID)
        nodes.push_back(maternal_id);

    if(paternal_id != UNKNOWN_ID)
        nodes.push_back(paternal_id);

    for(unsigned int i = 0; i < children.size(); ++i)
        nodes.push_back(children[i]->internal_id);

    for(unsigned int i = 0; i < mates.size(); ++i)
        nodes.push_back(mates[i]->internal_id);
        
    // childrens genotypes are independent of their siblings, but in
    // (highly?) consanguineous pedigrees you need to include them as
    // neighbours if they have already been peeled to get their half of 
    // the cutset that was peeled on to your parents
    
    // get mothers children
    // for each one where father is the same man
    // add to nodes if it is already peeled
    
    if(isfounder())
        return;
    
    Person* p = ped->get_by_index(maternal_id);
    for(unsigned i = 0; i < p->num_children(); ++i) {
        Person* tmp = p->get_child(i);
        if((tmp->get_paternalid() == paternal_id) and ps.is_peeled(tmp->internal_id)) {
            nodes.push_back(tmp->internal_id);
        }
    }
}
/*
bool Person::is_parent(unsigned int i) {
    return (i == maternal_id) or (i == paternal_id);
}
*/
void Person::get_cutset(PeelOperation& operation, PeelingState& state) {
    queue<unsigned int> q;
    vector<unsigned int> n;
    vector<int> visited(ped->num_members(), WHITE);
    unsigned int tmp, tmp2;
    Person* p;

    
    state.toggle_peel_operation(operation);
    
    visited[internal_id] = GREY;
    q.push(internal_id);

/*    
    // XXX this assumes that we are generating a simple peeling 
    // sequence on an arbitrary graph...
    for(unsigned int i = 0; i < ped->num_members(); ++i) {
        if(state.is_peeled(i)) {
            visited[i] = GREY;
            q.push(i);
            break;
        }
    }
*/
    while(not q.empty()) {
		tmp = q.front();
		q.pop();
		
        p = ped->get_by_index(tmp);
        p->neighbours(n, state);

        for(unsigned int i = 0; i < n.size(); ++i) {
            tmp2 = n[i];

            if(not state.is_peeled(tmp2)) {
                //if(tmp2 != internal_id) {
                    operation.add_cutnode(tmp2);
                //}
                continue;
            }

            if(visited[tmp2] == WHITE) {
                visited[tmp2] = GREY;
                q.push(tmp2);
            }
        }
        
        visited[tmp] = BLACK;
    }
    
    state.toggle_peel_operation(operation);
}

