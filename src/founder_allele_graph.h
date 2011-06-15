#ifndef LKG_FOUNDERALLELEGRAPH_H_
#define LKG_FOUNDERALLELEGRAPH_H_

#include <vector>

#include "genotype.h"


class Pedigree;
class DescentGraph;
class GeneticMap;

struct adj_node {
    int id;
    enum unphased_genotype label;
};

class FounderAlleleGraph {
	
    GeneticMap* map;
    Pedigree* ped;
    int num_founder_alleles;
    
    int* num_neighbours;
    adj_node** adj_matrix;
    
    int* best_descentstate;
    
    void _init();
    void _copy(const FounderAlleleGraph& rhs);
    void _kill();
    int _get_num_neighbours(int node) const;
    bool _add(int mat_fa, int pat_fa, enum unphased_genotype g);
    bool _check_legality(int node, int node_assignment, vector<int>& assignments);
    void _assign_and_recurse(int *component, int component_size, unsigned locus, 
                             int current_index, vector<int>& assignment, double *prob, double *best);
    double _enumerate_component(int *component, int component_size, unsigned locus);
    
 public :
	FounderAlleleGraph(GeneticMap* g, Pedigree* p);
	FounderAlleleGraph(const FounderAlleleGraph& fag);    
	~FounderAlleleGraph();
	FounderAlleleGraph& operator=(const FounderAlleleGraph& fag);
    
    void reset();    
	bool populate(DescentGraph& d, unsigned locus);
	bool likelihood(double* probs, unsigned locus);
    double descentstate_likelihood(unsigned locus);
    int get_founderallele_assignment(int fa);
    
    void print() const;
    void print_ds();
};

#endif

