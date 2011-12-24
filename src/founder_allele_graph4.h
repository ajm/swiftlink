#ifndef LKG_FOUNDERALLELEGRAPH4_H_
#define LKG_FOUNDERALLELEGRAPH4_H_

#include <vector>
#include "genotype.h"
#include "trait.h"
#include "genetic_map.h"


class Pedigree;
class DescentGraph;

#define DEFAULT_COMPONENT -1

class FounderAlleleGraph4 {
	
    Pedigree* ped;
    GeneticMap* map;
    unsigned int locus;
    unsigned int founder_alleles;
    double major_freq;
    double minor_freq;
    
    vector<int> edge_list;
    vector<int> group_membership;
    vector<int> group_fixed; // -1,0 or 1 - ie: -1 unfixed, 0,1 fixed to that index
    vector<vector<enum unphased_genotype> > allele_assignment;
    
    vector<int>* sequence;
    
    bool legal(enum unphased_genotype obs, enum unphased_genotype a1, enum unphased_genotype a2);
    enum unphased_genotype get_other_allele(enum unphased_genotype obs, enum unphased_genotype a1);
    void combine_components(int component1, int component2, bool flip);
    void propagate_fa_update(Person* p, int old_fa, int new_fa);
    double get_freq(enum unphased_genotype g);
    
 public :
	FounderAlleleGraph4(Pedigree* p, GeneticMap* g, int locus) :
        ped(p),
        map(g),
        locus(locus),
        founder_alleles(p->num_founders() * 2),
        major_freq(g->get_major(locus)),
        minor_freq(g->get_minor(locus)),
        edge_list(p->num_members() * 2, 0),
        group_membership(founder_alleles, DEFAULT_COMPONENT),
        group_fixed(founder_alleles, -1),
        allele_assignment(2, vector<enum unphased_genotype>(founder_alleles, UNTYPED)),
        sequence(NULL) {}
    
	FounderAlleleGraph4(const FounderAlleleGraph4& f) :
        ped(f.ped),
        map(f.map),
        locus(f.locus),
        founder_alleles(f.founder_alleles),
        major_freq(f.major_freq),
        minor_freq(f.minor_freq),
        edge_list(f.edge_list),
        group_membership(f.group_membership),
        group_fixed(f.group_fixed),
        allele_assignment(f.allele_assignment),
        sequence(f.sequence) {}
    
	~FounderAlleleGraph4() {}
    
	FounderAlleleGraph4& operator=(const FounderAlleleGraph4& rhs) {
        if(&rhs != this) {
            map = rhs.map;
            ped = rhs.ped;
            locus = rhs.locus;
            founder_alleles = rhs.founder_alleles;
            major_freq = rhs.major_freq;
            minor_freq = rhs.minor_freq;
            
            edge_list = rhs.edge_list;
            group_membership = rhs.group_membership;
            group_fixed = rhs.group_fixed;
            allele_assignment = rhs.allele_assignment;
            
            sequence = rhs.sequence;
        }
        
        return *this;
    }
    
    void set_sequence(vector<int>* seq) { sequence = seq; }
    
    double init_likelihood(DescentGraph& d, int locus);
    double update_likelihood(unsigned int personid, enum parentage p); // need DescentGraph as well?
    
    
    string debug_string();
};

#endif

