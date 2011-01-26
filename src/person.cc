using namespace std;

#include <cstdio>
#include <vector>
#include <queue>

#include "lkg.h"
#include "person.h"
#include "pedigree.h"
#include "genotype.h"
#include "peeling.h"


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
    return ps.is_peeled(this->maternal_id) and ps.is_peeled(this->paternal_id);
}

bool Person::ripe_above(PeelingState& ps) {
    return parents_peeled(ps);
}

bool Person::ripe_below(PeelingState& ps) {
    return offspring_peeled(ps);
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

bool Person::ripe(PeelingState& ps) {
    return ripe_all(ps) or \
            ripe_across(ps) or \
            ripe_below(ps) or \
            ripe_above_partners(ps);
}

bool Person::peel_operation(PeelOperation* p, PeelingState& state) {
    p->set_pivot(internal_id);
    
    if(ripe(state)) {
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

void Person::get_cutset(PeelOperation* operation, PeelingState& state) {
    queue<unsigned int> q;
    vector<unsigned int> n;
    unsigned int visited[ped->num_members()];
    unsigned int tmp, tmp2;
    Person* p;

    for(unsigned int i = 0; i < ped->num_members(); ++i) {
		visited[i] = WHITE;
    }
    
    visited[internal_id] = GREY;
    q.push(internal_id);

    while(not q.empty()) {
		tmp = q.front();
		q.pop();
		
        p = ped->get_by_index(tmp);
        p->neighbours(n);

        for(unsigned int i = 0; i < n.size(); ++i) {
            tmp2 = n[i];

            if(not state.is_peeled(tmp2)) {
                if(internal_id != tmp2) {
                    operation->add_cutnode(tmp2);
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
}

