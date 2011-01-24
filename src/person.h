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

class Person {
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
	
	vector<Person*> children;
	vector<Person*> mates;
	
	Pedigree* ped;
	
    // for peeling code
    //double disease_probability[4];

	bool _is_unknown(const string& s) const { return s == "0"; }
	string gender_str() const ;
	string affection_str() const ;

 public :
	Person(const string name, const string father_name, const string mother_name, 
			enum sex s, enum affection a, Pedigree* pedigree) :
			id(name),
			mother(mother_name),
			father(father_name),		
			gender(s),
			affection(a),			
			internal_id(UNKNOWN_ID),
			maternal_id(UNKNOWN_ID),
			paternal_id(UNKNOWN_ID),
			typed(true),
			ped(pedigree) {}
	
	~Person() {
		//fprintf(stderr, "Person %s destructor\n", id.c_str());
	}

	/* getters */
	string& get_id() { return this->id; }
	string& get_mother() { return this->mother; }
	string& get_father() { return this->father; }
	
	unsigned int get_internalid() const { return this->internal_id; }
	unsigned int get_maternalid() const { return this->maternal_id; }
	unsigned int get_paternalid() const { return this->paternal_id; }

	enum sex get_sex() const { return this->gender; }
	enum affection get_affection() const { return this->affection; }

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
	void set_internalid(unsigned int id) { this->internal_id = id; }
	void set_maternalid(unsigned int id) { this->maternal_id = id; }
	void set_paternalid(unsigned int id) { this->paternal_id = id; }
	void set_untyped() { this->typed = false; }

	void add_genotype(unphased_genotype_t g) {
		genotypes.push_back(g);
	}

//	void person_set_diseaseprob(person_t *p, config_t *c);

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
	
	bool mendelian_errors() const;
	void fill_in_relationships();
	void print() const;

	// so I can sort, I don't care for a specific ordering, 
	// I just want all the founders first
	bool operator<(const Person& p) const {
		return isfounder() and not p.isfounder();
	}

    bool peel_operation(PeelOperation* p);
};

#endif

