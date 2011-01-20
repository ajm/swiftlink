using namespace std;

#include <cstdio>
#include <string>
#include <algorithm>
#include <queue>

#include "lkg.h"
#include "genotype.h"
#include "pedigree.h"
#include "person.h"

bool Pedigree::sanity_check() {
    int components;

	if( not _same_number_of_markers() ) {
		return false;
	}
	
    if( _parental_relationship_errors() ) {
		//printf("parent relationship errors\n");
        return false;
	}
	
    _reorder_family();
	_rename_parent_ids();
	_fill_in_relationships();
	
    if( (components = _count_components()) != 1 ) {
        fprintf(stderr, "error: %s, family %s is actually composed of %d \
distinct families\n", __func__, id.c_str(), components);
        return false;
    }
	
    // ideally this should say what marker and what column, individuals
    if( _mendelian_errors() ) {
		//printf("mendelian errors\n");
        return false;
    }
	
	return true;
}

bool Pedigree::_same_number_of_markers() {
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
    for(uint i = 0; i < members.size(); ++i) {
        p = &members[i];
		
		id = p->get_id();
		
        // mother
        mother = p->get_mother();
		
        if( not p->mother_unknown() ) {
			if( (tmp = get_by_name(mother)) == NULL ) {
                fprintf(stderr, "error: %s, mother of %s does not exist\n", 
						__func__, id.c_str());
                error = true;
            }
            else if( not tmp->isfemale() ) {
                fprintf(stderr, "error: %s, mother of %s is not female\n", 
						__func__, id.c_str());
                error = true;
            }
        }
		
        // father
        father = p->get_father();
        
		if( not p->father_unknown() ) {
            if( (tmp = get_by_name(father)) == NULL ) {
                fprintf(stderr, "error: %s, father of %s does not exist\n", 
						__func__, id.c_str());
                error = true;
            }
            else if( not tmp->ismale() ) { 
                fprintf(stderr, "error: %s, father of %s is not male\n", 
						__func__, id.c_str());
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
    for(uint i = 0; i < members.size(); ++i) {
        members[i].set_internalid(i);
    }
}

void Pedigree::_rename_parent_ids() {
	Person* p;
	Person* tmp;
	string id, mother, father;
	
	// each person, look at maternal and paternal ids and search
	// the family for that persons id to find out their internal id
	for(uint i = 0; i < members.size(); ++i) {
		p = &members[i];
		
		id	   = p->get_id();
		mother = p->get_mother();
		father = p->get_father();
		
		if( p->mother_unknown() )
			p->set_maternalid(UNKNOWN_PARENT);
		
		if( p->father_unknown() )
			p->set_paternalid(UNKNOWN_PARENT);
		
		for(uint j = 0; j < members.size(); ++j) {
			tmp = &members[j];
            
			if(tmp->get_id() == mother)
				p->set_maternalid(j);

			if(tmp->get_id() == father)
				p->set_paternalid(j);
		}
	}
}

bool Pedigree::_mendelian_errors() const {
	for(uint i = 0; i < members.size(); ++i) {
		if(members[i].mendelian_errors())
			return true;
	}
    return false;
}

void Pedigree::_fill_in_relationships() {
	for(uint i = 0; i < members.size(); ++i) {
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
	int visited[members.size()];
	Person* p;

	for(uint i = 0; i < members.size(); ++i)
		visited[i] = WHITE;

	do {
		// find a starting point
		for(uint i = 0; i < members.size(); ++i) {
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

			for(uint i = 0; i < p->num_children(); ++i) {
				tmp2 = p->get_child(i)->get_internalid();
				if(visited[tmp2] == WHITE) {
					visited[tmp2] = GREY;
					q.push(tmp2);
				}
			}

			tmp2 = p->get_maternalid();
			if((tmp2 != UNKNOWN_PARENT) && (visited[tmp2] == WHITE)) {
                visited[tmp2] = GREY;
                q.push(tmp2);
            }

			tmp2 = p->get_paternalid();
			if((tmp2 != UNKNOWN_PARENT) && (visited[tmp2] == WHITE)) {
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

void Pedigree::print() const {
	printf("Pedigree: %s\n", id.c_str());
	for(uint i = 0; i < members.size(); ++i) {
		members[i].print();
		printf("\n");
	}
}

