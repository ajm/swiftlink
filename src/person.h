#ifndef LKG_PERSON_H_
#define LKG_PERSON_H_

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "genotype.h"

const unsigned int UNKNOWN_PARENT = ~0u;
const unsigned int UNKNOWN_ID = UNKNOWN_PARENT;

enum sex {
	UNSEXED,
    MALE,
    FEMALE
};

enum affection {
	UNKNOWN_AFFECTION,
    UNAFFECTED,
    AFFECTED
};

class Pedigree;
class PeelOperation;
class PeelingState;
class DiseaseModel;

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
	vector<unphased_genotype_t> genotypes;
	
    // filled in by the pedigree class once
    // the entire pedigree file has been read
	vector<Person*> children;
	vector<Person*> mates;
    
    // allows me to get family members as a 
    // call to a person object
	Pedigree* ped;
	
    // for peeling code, information comes from the
    // disease model objects
    double disease_prob[4];


    // private stuff
	bool _is_unknown(const string& s) const { 
        return s == "0"; 
    }
	string gender_str() const;
	string affection_str() const;


 public :
	Person(const string name, const string father_name, const string mother_name, 
			enum sex s, enum affection a, Pedigree* pedigree/*, DiseaseModel& dm*/) :
			id(name),
			mother(mother_name),
			father(father_name),		
			gender(s),
			affection(a),			
			internal_id(UNKNOWN_ID),
			maternal_id(UNKNOWN_ID),
			paternal_id(UNKNOWN_ID),
			typed(true),
			ped(pedigree) {
        
//        init_probs(dm);
    }
	
	~Person() {}

	/* getters */
	string& get_id() { return id; }
	string& get_mother() { return mother; }
	string& get_father() { return father; }
	
	unsigned int get_internalid() const { return internal_id; }
	unsigned int get_maternalid() const { return maternal_id; }
	unsigned int get_paternalid() const { return paternal_id; }

	enum sex get_sex() const { return gender; }
	enum affection get_affection() const { return affection; }

	unphased_genotype_t get_genotype(unsigned int i) const {
		if(not istyped())
			return UNTYPED;
		
#ifndef _FEELING_LUCKY_
		if(genotypes.size() < i)
			abort();
#endif
		
		return genotypes[i];
	}
	
	unsigned int num_markers() { return genotypes.size(); }
	unsigned int num_children() { return children.size(); }
	unsigned int num_mates() { return mates.size(); }
	Person* get_child(unsigned int i) const { return children[i]; }
	Person* get_spouse(unsigned int i) const { return mates[i]; }
	Person* get_mate(unsigned int i) const { return mates[i]; }
	unphased_genotype_t get_marker(unsigned i) const { return genotypes[i]; }

	/* setters */
	void set_internalid(unsigned int id) { internal_id = id; }
	void set_maternalid(unsigned int id) { maternal_id = id; }
	void set_paternalid(unsigned int id) { paternal_id = id; }
	void set_untyped() { typed = false; }

	void add_genotype(unphased_genotype_t g) {
		genotypes.push_back(g);
	}

	/* tests */
	bool ismale() const { return this->gender == MALE; }
	bool isfemale() const { return this->gender == FEMALE; }
	bool isunsexed() const { return this->gender == UNSEXED; }

	bool isaffected() const { return this->affection == AFFECTED; }
	bool istyped() const { return this->typed; }

	bool isfounder() const { return mother_unknown() and father_unknown(); }
	bool isleaf() const { return children.size() == 0; }

	bool mother_unknown() const { return _is_unknown(mother); }
	bool father_unknown() const { return _is_unknown(father); }
	
    // pedigree construction / validation
	bool mendelian_errors() const;
	void fill_in_relationships();
	void print() const;

	// so I can sort, I don't care for a specific (strong) ordering, 
	// I just want all the founders first
	bool operator<(const Person& p) const {
		return isfounder() and not p.isfounder();
	}

    // peeling
    // XXX many of these can be private
    bool contains_unpeeled(vector<Person*>& v, PeelingState& ps);
    unsigned int count_unpeeled(vector<Person*>& v, PeelingState& ps);
    bool offspring_peeled(PeelingState& ps);
    bool partners_peeled(PeelingState& ps);
    bool parents_peeled(PeelingState& ps);
    bool ripe_above(PeelingState& ps);
    bool ripe_below(PeelingState& ps);
    bool ripe_across(PeelingState& ps);
    bool ripe_all(PeelingState& ps);
    bool ripe_above_partners(PeelingState& ps);
    bool ripe(PeelingState& ps);
    bool peel_operation(PeelOperation* p, PeelingState& state);
    void neighbours(vector<unsigned int>& nodes);
    void get_cutset(PeelOperation* operation, PeelingState& state);

    // XXX should be private and disease model passed as an argument
    // to constructor
    void init_probs(DiseaseModel& dm);
};

#endif

