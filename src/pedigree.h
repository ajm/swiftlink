#ifndef LKG_PEDIGREE_H_
#define LKG_PEDIGREE_H_

using namespace std;

#include <cstdio>
#include <string>
#include <vector>

#include "types.h"
#include "person.h"


class Pedigree {

	string id;
	vector<Person> members;
	unsigned int number_of_founders;
	unsigned int number_of_leaves;
	
	
	bool _mendelian_errors() const;
	bool _parental_relationship_errors();
	bool _same_number_of_markers() const;
    void _reorder_family();
	void _rename_parent_ids();
	void _fill_in_relationships();
    void _set_typed_flags();
	void _count_founders();
    void _count_leaves();
	int _count_components();
	int _person_compare(const Person& a, const Person& b) const;
    
 public:
	Pedigree(const string id) : 
	    id(id), 
        members(), 
        number_of_founders(0), 
        number_of_leaves(0) {}
    
    Pedigree(const Pedigree& rhs) :
        id(rhs.id),
        members(rhs.members),
        number_of_founders(rhs.number_of_founders),
        number_of_leaves(rhs.number_of_leaves) {}
    
	~Pedigree() {}

    Pedigree& operator=(const Pedigree& rhs) {
        
        if(&rhs != this) {
            id = rhs.id;
            members.clear();
            members = rhs.members;
            number_of_founders = rhs.number_of_founders;
            number_of_leaves = rhs.number_of_leaves;
        }
        
        return *this;
    }

	// interrogate	
	string get_id() const { 
        return id;
    }
    
	unsigned int num_members() const { 
        return members.size();
    }
    
   	unsigned int num_markers() const {
		return members[0].num_markers(); 
	}
    
	unsigned int num_founders() const {
		return number_of_founders;
	}
    
	unsigned int num_leaves() const {
		return number_of_leaves;
	}
    
	//Person* get_by_index(int i);
	Person* get_by_name(const string& id);
	
	inline Person* get_by_index(int i) {
        return &members[i];
    }

	
	// manipulate
	bool add(Person& p);
	bool exists(const string& id);

	bool sanity_check();
	string debug_string();
};

#endif

