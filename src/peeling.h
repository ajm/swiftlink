#ifndef LKG_GENOTYPE_H_
#define LKG_GENOTYPE_H_

using namespace std;

#include <vector>

#include "pedigree.h"

class PeelOperation {
    
    unsigned int pivot;
    vector<unsigned int>* cutset;
    
  public :
    PeelOperation(unsigned int pivot_node, vector<unsigned int>* cutset_nodes) 
        : pivot(pivot_node), cutset(cutset_nodes) {}
    
    unsigned int get_pivot() { 
        return pivot;
    }
    
    unsigned int get_cutset_size() { 
        return cutset.size();
    }

    void set_pivot(unsigned int p) {
        pivot = p;
    }

    void set_cutset(vector<unsigned int>* c) {
        cutset = c;
    }
};

// XXX i think I would need one that could peel :
// i) over genotypes at loci (to create L-sampler)
// ii) over descent graphs in-between loci (to calculate the LOD score)
//
// notes: peeling sequences might not be the same, because the number of
// things to peel over is different, ie: 2 alleles make 3 genotypes, but
// they make 4 possible descent state settings
// (if my understanding of this is correct, genotypes are unphased, but
// the descent graphs require phase)

class PedigreePeeler {

    Pedigree& ped;
    vector<PeelOperation> peel;
    bool* peeled;

    PeelOperation get_minimum_cutset(vector<PeelOperation*> peels);
    vector<PeelOperation> get_possible_peels(unsigned int* unpeeled);

  public :
    PedigreePeeler(Pedigree& p) : ped(p) {
        peeled = new bool[ped.size()];
    }

    ~PedigreePeeler() {
        delete[] peeled;
    }

    bool build_peel_order();
    void print_peel_order();
    bool peel(double *likelihood);
};

#endif

