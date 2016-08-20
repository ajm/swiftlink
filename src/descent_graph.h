#ifndef LKG_DESCENTGRAPH_H_
#define LKG_DESCENTGRAPH_H_

using namespace std;

#include <limits>
#include <string>

#include "types.h"
#include "logarithms.h"
#include "genetic_map.h"

class Pedigree;

class DescentGraph {

	int* data;
	Pedigree* ped;
	GeneticMap* map;
	double prob;
    double marker_transmission; // cache for transmission prob
    int graph_size; 			// size of descent graph at one loci, 
								// for indexing data
	int recombinations;
    bool sex_linked;
	
	vector<int> seq;
	
    void _invalidate_paternal_x();
	double _transmission_prob();
	double _recombination_prob();
    double _best_prior_prob();
	double _sum_prior_prob();
    //int _offset(unsigned person_id, unsigned locus, enum parentage p) const ;
	inline int _offset(unsigned person_id, unsigned locus, enum parentage p) const {
	    return (graph_size * locus) + (person_id * 2) + p;
    }
	int _founder_allele(unsigned person_id, enum parentage p) const;
	void find_founderallelegraph_ordering();

 public :
	DescentGraph(Pedigree* ped, GeneticMap* map, bool sex_linked);
	DescentGraph(const DescentGraph& d);
    ~DescentGraph();
	DescentGraph& operator=(const DescentGraph& d);
    
    bool operator<(const DescentGraph& a) const {
        return prob < a.prob;
    }
	
	//void copy_from(DescentGraph& d, unsigned start, unsigned end);
    
    //int get(unsigned person_id, unsigned locus, enum parentage p) const ;
    inline int get(unsigned person_id, unsigned locus, enum parentage p) const {
        return data[_offset(person_id, locus, p)];
    }
    void set(unsigned person_id, unsigned locus, enum parentage p, int value);
    void flip(unsigned person_id, unsigned locus, enum parentage p);
    
    int get_bit(unsigned i) const ;
    void set_bit(unsigned i, int b);
    void flip_bit(unsigned i);
	
    void copy_locus(unsigned src, unsigned dst);
    void copy_locus(DescentGraph& d, unsigned src, unsigned dst);

    int get_founderallele(unsigned person_id, unsigned loci, enum parentage p) const;
    
    bool random_descentgraph();
    bool illegal() const { return prob == LOG_ILLEGAL; }
	int num_recombinations() const { return recombinations; }
    double get_prob() const { return prob; }
	double get_marker_transmission() const { return marker_transmission; }
	double get_recombination_prob(unsigned int locus, bool count_crossovers);
    double get_haplotype_likelihood();
    double get_likelihood();

    double get_likelihood2(GeneticMap* m) {
        GeneticMap* tmp = map;
        map = m;
        double l = get_likelihood();
        map = tmp;
        return l;
    }
    
    string debug_string();
    
    int* get_internal_ptr() { return data; }
    size_t get_internal_size() { 
        return sizeof(int) * graph_size * map->num_markers();
    }
};

#endif

