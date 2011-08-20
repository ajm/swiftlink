#ifndef LKG_PERSON_H_
#define LKG_PERSON_H_

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "types.h"


class Pedigree;
class DiseaseModel;
class PeelOperation;
class PeelingState;

class Person {
    
    // info directly from pedigree file, 
    // essentially verbatim
	string id;
	string mother;
	string father;
	
	enum sex gender;
	enum affection affection;

	unsigned int internal_id;
	unsigned int maternal_id;
	unsigned int paternal_id;
	
	bool typed;
    
    // allows me to get family members as a 
    // call to a person object
	Pedigree* ped;
	
    // for peeling code, information comes from the
    // disease model objects
    double disease_prob[4];
    
    // filled in by the pedigree class once
    // the entire pedigree file has been read
	vector<enum unphased_genotype> genotypes;
	vector<Person*> children;
	vector<Person*> mates;


    
    // private stuff
	void _init_probs(const DiseaseModel& dm);
	bool _is_unknown(const string& s) const { 
        return s == "0"; 
    }

    // private, peeling related
    unsigned count_unpeeled(vector<Person*>& v, PeelingState& ps);
    unsigned get_unpeeled_mate(PeelingState& ps);
    bool offspring_peeled(PeelingState& ps);
    bool founder_mates_peeled(PeelingState& ps);
    bool partners_peeled(PeelingState& ps);
    bool parents_peeled(PeelingState& ps);
    bool one_parent_peeled(PeelingState& ps);
    bool ripe_above_singular_mate(PeelingState& ps);
    bool ripe_above(PeelingState& ps);
    bool ripe_above_at_least_one_parent(PeelingState& ps);
    
    bool ripe_to_peel_across(PeelingState& ps);
    bool ripe_to_peel_final(PeelingState& ps);
    bool ripe_to_peel_down(PeelingState& ps);
    bool ripe_to_peel_up(PeelingState& ps);
    
    void neighbours(vector<unsigned int>& nodes, PeelingState& ps);
    void get_cutset(PeelOperation& operation, PeelingState& state);
    void add_mate(Person* p);
    
    
    
 public :
	Person(const string name, const string father_name, const string mother_name, 
			enum sex s, enum affection a, Pedigree* pedigree, const DiseaseModel& dm);
	~Person() {}
	Person(const Person& rhs);
	Person& operator=(const Person& p);
    
    
    /* getters */
	string get_id() const { return id; }
	string get_mother() const { return mother; }
	string get_father() const { return father; }
	unsigned int get_internalid() const { return internal_id; }
	unsigned int get_maternalid() const { return maternal_id; }
	unsigned int get_paternalid() const { return paternal_id; }
	enum sex get_sex() const { return gender; }
	enum affection get_affection() const { return affection; }

	enum unphased_genotype get_genotype(unsigned int i) const {
		if(not istyped()) {
			return UNTYPED;
        }
        
		return genotypes[i];
	}
	enum unphased_genotype get_marker(unsigned i) const { return genotypes[i]; } // <--- redundancy?
    
	unsigned int num_markers() const { return genotypes.size(); }
	unsigned int num_children() const { return children.size(); }
	unsigned int num_mates() const { return mates.size(); }
	Person* get_child(unsigned int i) const { return children[i]; }
	//Person* get_spouse(unsigned int i) const { return mates[i]; }
	//Person* get_mate(unsigned int i) const { return mates[i]; }
	

	/* setters */
	void set_internalid(unsigned int id) { internal_id = id; }
	void set_maternalid(unsigned int id) { maternal_id = id; }
	void set_paternalid(unsigned int id) { paternal_id = id; }
    
	void add_genotype(enum unphased_genotype g) {
		genotypes.push_back(g);
	}

	/* tests */
	bool ismale() const { return gender == MALE; }
	bool isfemale() const { return gender == FEMALE; }
	bool isunsexed() const { return gender == UNSEXED; }

	bool isaffected() const { return affection == AFFECTED; }
	bool istyped() const { return typed; }

	bool isfounder_str() const { return mother_unknown() and father_unknown(); }
    bool isfounder() const /*{ return isfounder_str(); }*/ { return (maternal_id == UNKNOWN_PARENT) and (paternal_id == UNKNOWN_PARENT); }
	bool isleaf() const { return children.size() == 0; }

	bool mother_unknown() const { return _is_unknown(mother); }
	bool father_unknown() const { return _is_unknown(father); }
	
    // pedigree construction / validation
	bool mendelian_errors() const ;
	void fill_in_relationships();
    void set_typed();
    
	// so I can sort, I don't care for a specific (strong) ordering, 
	// I just want all the founders first
	bool operator<(const Person& p) const {
		return isfounder_str() and not p.isfounder_str();
	}
    
    double get_disease_prob(enum phased_trait pt) { return disease_prob[pt]; }
    //bool is_parent(unsigned int i);
    inline bool is_parent(unsigned int i) const {
        return (i == maternal_id) or (i == paternal_id);
    }
    
    bool peel_operation(PeelOperation& p, PeelingState& state);
    
    string debug_string();
};

#endif
