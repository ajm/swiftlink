using namespace std;

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <queue>

#include "types.h"
#include "pedigree.h"
#include "person.h"


bool Pedigree::add(Person& p) {
    if(not exists(p.get_id())) {
        members.push_back(p);
		return true;
    }

    return false;
}

bool Pedigree::exists(const string& id) {
	return get_by_name(id) != NULL;
}

/*
Person* Pedigree::get_by_index(int i) {
    return &members[i];
}
*/

Person* Pedigree::get_by_name(const string& id) {
	for(unsigned int i = 0; i < members.size(); ++i) {
		if(members[i].get_id() == id) {
			return &members[i];
		}
	}
	return NULL;
}

void Pedigree::_count_founders() {
	number_of_founders = 0;
		
	for(unsigned int i = 0; i < members.size(); ++i) {
		if(members[i].isfounder()) {
			number_of_founders++;
		}
	}
}

void Pedigree::_count_leaves() {
	number_of_leaves = 0;
	
	for(unsigned int i = 0; i < members.size(); ++i) {
		if(members[i].isleaf()) {
			number_of_leaves++;
		}
	}
}

bool Pedigree::sanity_check() {
    int components;

	if(not _same_number_of_markers())
		return false;
	
    if(_parental_relationship_errors())
        return false;
	
    _reorder_family();
	_rename_parent_ids();
	_fill_in_relationships();
	
    if((components = _count_components()) != 1) {
        fprintf(stderr, "error: %s, family %s is actually composed of %d distinct families\n", 
            __func__, id.c_str(), components);
        return false;
    }
	
    if(_mendelian_errors()) {
        return false;
    }
    
    _count_founders();
    _count_leaves();
	
	return true;
}

bool Pedigree::_same_number_of_markers() const {
	if(int(members.size()) < 2)
		return true;

	unsigned int tmp = members[0].num_markers();

	for(int i = 1; i < int(members.size()); ++i) {
		if(tmp != members[i].num_markers())
			return false;
	}

	return true;
}

bool Pedigree::_parental_relationship_errors() {
	bool error;	
    Person *p, *tmp;
	string id;
	string mother; 
	string father;
	
	error = false;
	
    // check that parents exist
    // check that parents are the correct genders
    for(unsigned int i = 0; i < members.size(); ++i) {
        p = &members[i];
		
		id = p->get_id();
		
        // mother
        mother = p->get_mother();
		
        if( not p->mother_unknown() ) {
			if( (tmp = get_by_name(mother)) == NULL ) {
                fprintf(stderr, "error: %s, mother of %s does not exist\n", __func__, id.c_str());
                error = true;
            }
            else if( not tmp->isfemale() ) {
                fprintf(stderr, "error: %s, mother of %s is not female\n", __func__, id.c_str());
                error = true;
            }
        }
		
        // father
        father = p->get_father();
        
		if( not p->father_unknown() ) {
            if( (tmp = get_by_name(father)) == NULL ) {
                fprintf(stderr, "error: %s, father of %s does not exist\n", __func__, id.c_str());
                error = true;
            }
            else if( not tmp->ismale() ) { 
                fprintf(stderr, "error: %s, father of %s is not male\n", __func__, id.c_str());
                error = true;
            }
        }
    }
	
    return error;
}

void Pedigree::_reorder_family() {
    // reorder so founders are from 0->n
	sort(members.begin(), members.end());
    
	// rename so id is still index
    for(unsigned int i = 0; i < members.size(); ++i) {
        members[i].set_internalid(i);
    }
}

void Pedigree::_rename_parent_ids() {
	Person* p;
	Person* tmp;
	string id, mother, father;
	
	// each person, look at maternal and paternal ids and search
	// the family for that persons id to find out their internal id
	for(unsigned int i = 0; i < members.size(); ++i) {
		p = &members[i];
		
		id	   = p->get_id();
		mother = p->get_mother();
		father = p->get_father();
		
		if(p->mother_unknown())
			p->set_maternalid(UNKNOWN_PARENT);
		
		if(p->father_unknown())
			p->set_paternalid(UNKNOWN_PARENT);
		
		for(unsigned int j = 0; j < members.size(); ++j) {
			tmp = &members[j];
            
			if(tmp->get_id() == mother)
				p->set_maternalid(j);

			if(tmp->get_id() == father)
				p->set_paternalid(j);
		}
	}
}

bool Pedigree::_mendelian_errors() const {
	for(unsigned int i = 0; i < members.size(); ++i) {
		if(members[i].mendelian_errors())
			return true;
	}
    return false;
}

void Pedigree::_fill_in_relationships() {
	for(unsigned int i = 0; i < members.size(); ++i) {
		members[i].fill_in_relationships();
	}
}

// pedigree BFS to assess connectedness and return number of distinct graphs
// obviously should be 1 if it is a correct family
int Pedigree::_count_components() {
	int components = 0;
	int total = members.size();
	unsigned int tmp, tmp2;
	queue<unsigned int> q;
    Person* p;
    vector<int> visited(members.size(), WHITE);

	do {
		// find a starting point
		for(unsigned int i = 0; i < members.size(); ++i) {
			if(visited[i] == WHITE) {
				q.push(i);
				visited[i] = GREY;
				break;
			}
		}

		while(not q.empty()) {
			tmp = q.front();
			q.pop();
			
			p = &members[tmp];

			for(unsigned int i = 0; i < p->num_children(); ++i) {
				tmp2 = p->get_child(i)->get_internalid();
				if(visited[tmp2] == WHITE) {
					visited[tmp2] = GREY;
					q.push(tmp2);
				}
			}

			tmp2 = p->get_maternalid();
			if((tmp2 != UNKNOWN_PARENT) and (visited[tmp2] == WHITE)) {
                visited[tmp2] = GREY;
                q.push(tmp2);
            }

			tmp2 = p->get_paternalid();
			if((tmp2 != UNKNOWN_PARENT) and (visited[tmp2] == WHITE)) {
                visited[tmp2] = GREY;
                q.push(tmp2);
            }

			visited[tmp] = BLACK;
			--total;
		}

		++components;

	} while(total != 0);
    
    
	return components;
}

string Pedigree::debug_string() {
    stringstream ss;
    
    ss << "Pedigree: " << id << endl;
    
    for(unsigned i = 0; i < members.size(); ++i) {
        ss << members[i].debug_string() << endl;
    }
        
    return ss.str();
}

