#ifndef LKG_FOUNDERALLELEGRAPH3_H_
#define LKG_FOUNDERALLELEGRAPH3_H_

#include <vector>
#include "genotype.h"


class Pedigree;
class DescentGraph;
class GeneticMap;

struct AdjacentNode {
    int id;
    enum unphased_genotype label;
};

class FounderAlleleGraph3 {
	
    Pedigree* ped;
    GeneticMap* map;
    unsigned int founder_alleles;
    
    vector<unsigned int> num_neighbours;
    vector< vector<AdjacentNode> > adj_matrix;
    vector<int> assignment;
    vector<int>* sequence;
    
    int get_num_neighbours(int node) const {
        return num_neighbours[node];
    }
    
    bool correct_alleles_loop(enum unphased_genotype g, int allele1);
    bool correct_alleles(enum unphased_genotype g, int allele1);
    bool correct_alleles(enum unphased_genotype g, int allele1, int allele2);
    bool legal(vector<int>& component, vector<int>& q);
    void assign(int allele1, int allele2, enum unphased_genotype g);
    bool add(int mat, int pat, enum unphased_genotype g);
    bool populate(DescentGraph& d, int locus);
    double enumerate(int locus, vector<int>& component);
    double component_likelihood(int locus, vector<int>& component);
    double likelihood(int locus);
    
 public :
	FounderAlleleGraph3(Pedigree* p, GeneticMap* g) :
        ped(p),
        map(g),
        founder_alleles(p->num_founders() * 2),
        num_neighbours(founder_alleles, 0),
        adj_matrix(founder_alleles), 
        assignment(ped->num_members() * 2),
        sequence(NULL) {
        
        for(unsigned i = 0; i < founder_alleles; ++i) {
            adj_matrix[i].resize(founder_alleles);
        }
    }
    
	FounderAlleleGraph3(const FounderAlleleGraph3& f) :
        ped(f.ped),
        map(f.map),
        founder_alleles(f.founder_alleles),
        num_neighbours(f.num_neighbours),
        adj_matrix(f.adj_matrix), 
        assignment(f.assignment),
        sequence(f.sequence) {}
    
	~FounderAlleleGraph3() {}
    
	FounderAlleleGraph3& operator=(const FounderAlleleGraph3& rhs) {
        if(&rhs != this) {
            map = rhs.map;
            ped = rhs.ped;
            
            founder_alleles = rhs.founder_alleles;
            
            num_neighbours = rhs.num_neighbours;
            adj_matrix = rhs.adj_matrix;
            assignment = rhs.assignment;
            sequence = rhs.sequence;
        }
        
        return *this;
    }
    
    void set_sequence(vector<int>* seq) { sequence = seq; }
    
    double evaluate(DescentGraph& d, unsigned locus);
    
    string debug_string();
};

#endif

