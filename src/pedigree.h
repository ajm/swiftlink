#ifndef LKG_PEDIGREE_H_
#define LKG_PEDIGREE_H_

using namespace std;

#include <cstdio>
#include <string>
#include <vector>

#include "lkg.h"
#include "person.h"


class Pedigree {

	string id;
	vector<Person> members;
	unsigned int number_of_founders;
	unsigned int number_of_leaves;
	
	
	bool _mendelian_errors() const;
	bool _parental_relationship_errors();
	bool _same_number_of_markers();
	int _count_components();
	
	unsigned int _count_founders() {
		number_of_founders = 0;
		
		for(uint i = 0; i < members.size(); ++i) {
			if(members[i].isfounder()) {
				number_of_founders++;
			}
		}

		return number_of_founders;
	}

	unsigned int _count_leaves() {
		number_of_leaves = 0;
		
		for(uint i = 0; i < members.size(); ++i) {
			if(members[i].isleaf()) {
				number_of_leaves++;
			}
		}

		return number_of_leaves;
	}
	
	int _person_compare(const Person& a, const Person& b) const;
	void _reorder_family();
	void _rename_parent_ids();
	void _fill_in_relationships();
	
 public:
	Pedigree(const string id) : id(id), number_of_founders(0), number_of_leaves(0) {}
	~Pedigree() {}

	// interrogate	
	string& get_id() { return this->id; }
	unsigned int num_members() const { return members.size(); }
	unsigned int num_founders() { 
		return number_of_founders != 0 ? 
			number_of_founders : 
			_count_founders();
	}
	unsigned int num_leaves() { 
		return number_of_leaves != 0 ? 
			number_of_leaves : 
			_count_leaves();
	}
	unsigned int num_markers() { 
		return members[0].num_markers(); 
	}
		
	Person* get_by_index(int i) { return &members[i]; }
	Person* get_by_name(const string& id) {
		for(uint i = 0; i < members.size(); ++i) {
			if(members[i].get_id() == id) {
				return &members[i];
			}
		}
		return NULL;
	}
	
	// manipulate
	bool add(Person* p) {
		for(uint i = 0; i < members.size(); ++i) {
			if(members[i].get_id() == p->get_id())
				return false;
		}
        
		members.push_back(*p);
		return true;
	}
	bool exists(const string& id) {
		return get_by_name(id) != NULL;
	}
		
	bool sanity_check();
	void print() const;
};

#endif

