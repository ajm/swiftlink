#ifndef LKG_PERSON_H_
#define LKG_PERSON_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "types.h"
#include "genotype.h"

using namespace std;


class Pedigree;
class DiseaseModel;
class PeelOperation;
class PeelingState;
class GeneticMap;

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
    DiseaseModel* dm;

    // filled in by the pedigree class once
    // the entire pedigree file has been read
	vector<enum unphased_genotype> genotypes;
    vector<vector<double> > genotypes_prob;
	vector<Person*> children;
	vector<Person*> mates;

    // private stuff
	void _init_probs();
	bool _is_unknown(const string& s) const { 
        return s == "0"; 
    }
    
    void add_mate(Person* p);
    unsigned count_unpeeled(vector<Person*>& v, PeelingState& ps);
    
    double get_genotype_probability(enum unphased_genotype g, enum phased_trait pt, double marker_prob);
    
 public :
	Person(const string name, const string father_name, const string mother_name, 
			enum sex s, enum affection a, Pedigree* pedigree, DiseaseModel& dm);
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
	
	unsigned int get_parentid(enum parentage p) const {
	    switch(p) {
	        case MATERNAL:
	            return maternal_id;
	        case PATERNAL:
	            return paternal_id;
            case NONE :
            default :
                break;
	    }
        abort();
	}

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
	Person* get_mate(unsigned int i) const { return mates[i]; }
	
	bool is_offspring(unsigned int node);

	/* setters */
	void set_internalid(unsigned int id) { internal_id = id; }
	void set_maternalid(unsigned int id) { maternal_id = id; }
	void set_paternalid(unsigned int id) { paternal_id = id; }
    
	void add_genotype(enum unphased_genotype g) {
	    if(g != UNTYPED) {
	        typed = true;
	    }
	    
		genotypes.push_back(g);
    }

    void clear_genotypes() {
        genotypes.clear();
    }
	
	bool legal_genotype(unsigned locus, enum phased_trait g) {
	
	    //fprintf(stderr, "id = %d, locus = %d, g = %d, typed = %s\n", internal_id, locus, g, typed ? "true" : "false");
	    
	    if(not typed)
	        return true;
	    
	    switch(genotypes[locus]) {
	        case UNTYPED :
	            return true;
	        case HETERO :
	            return (g == TRAIT_AU) or (g == TRAIT_UA);
	        case HOMOZ_A :
	            return g == TRAIT_UU;
	        case HOMOZ_B :
	            return g == TRAIT_AA;
	    }
	    
	    abort();
	}

	/* tests */
	bool ismale() const { return gender == MALE; }
	bool isfemale() const { return gender == FEMALE; }
	bool isunsexed() const { return gender == UNSEXED; }

	bool isaffected() const { return affection == AFFECTED; }
	void make_unknown_affection() { 
        affection = UNKNOWN_AFFECTION; 
        _init_probs();
    }
    bool istyped() const { return typed; }

	bool isfounder_str() const { return mother_unknown() and father_unknown(); }
    bool isfounder() const /*{ return isfounder_str(); }*/ { return (maternal_id == UNKNOWN_PARENT) and (paternal_id == UNKNOWN_PARENT); }
	bool isleaf() const { return children.size() == 0; }

	bool mother_unknown() const { return _is_unknown(mother); }
	bool father_unknown() const { return _is_unknown(father); }
	
    // pedigree construction / validation
	bool mendelian_errors() const ;
	void fill_in_relationships();
    
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
    
    bool partners_peeled(PeelingState& ps);
    bool offspring_peeled(PeelingState& ps);
    
    bool safe_to_ignore_meiosis(enum parentage p);
    
    string debug_string();

    void populate_trait_prob_cache(GeneticMap& map);
    double get_trait_probability(unsigned int locus, enum phased_trait pt) {
        return genotypes_prob[locus][pt];
    }

    // experimental for the ELOD code
    // assumes the genotype is already there
    void copy_disease_probs(int locus) {
        for(int i = 0; i < 4; ++i) {
            genotypes_prob[locus][i] = disease_prob[i];
        }
    }
};

#endif

