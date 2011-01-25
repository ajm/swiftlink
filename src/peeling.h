#ifndef LKG_PEELING_H_
#define LKG_PEELING_H_

using namespace std;

#include <vector>
#include <cstdio>

#include "pedigree.h"

class PeelingState {
    vector<bool> peeled;

    PeelingState(Pedigree* p) peeled(p->num_members(), false) {}

    bool is_peeled(unsigned int i) {
        return peeled[i];
    }

    void set_peeled(unsigned int i) {
        peeled[i] = true;
    }
};

class PeelOperation {
    
    unsigned int pivot;
    vector<unsigned int> cutset;
    
  public :
    PeelOperation() pivot(-1), cutset(0) {}
    
    PeelOperation(unsigned int pivot_node, vector<unsigned int>* cutset_nodes) 
        : pivot(pivot_node), cutset(cutset_nodes) {}

    ~PeelOperation() {
        delete cutset;
    }
    
    unsigned int get_pivot() const { 
        return pivot;
    }
    
    unsigned int get_cutset_size() const { 
        return cutset.size();
    }

    void set_pivot(unsigned int p) {
        pivot = p;
    }

    void add_cutnode(unsigned int c) {
        cutset.push_back(c);
    }

    void print() const {
        int tmp = cutset.size();

        printf("%d:");
        for(unsigned int i = 0; i < tmp; ++i) {
            printf("%d");
            if(i != (tmp-1)) {
                putchar(',');
            }
        }
        putchar('\n');
    }

    bool operator<(const PeelOperation& p) const {
		return cutset.size() < p.cutset.size();
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

    Pedigree* ped;
    vector<PeelOperation> peel;
    PeelingState state;

    PeelOperation get_minimum_cutset(vector<PeelOperation*> peels);
    vector<PeelOperation> get_possible_peels(unsigned int* unpeeled);

  public :
    PedigreePeeler(Pedigree* p) 
        : ped(p), ps(p) {}
    ~PedigreePeeler() {}

    PeelOperation get_random_operation(
            vector<PeelOperation>::iterator start,
            vector<PeelOperation>::iterator end
        );
    
    PeelOperation get_best_operation_heuristic(
            vector<PeelOperation>::iterator start,
            vector<PeelOperation>::iterator end
        );
    
    PeelOperation get_best_operation(
            vector<PeelOperation>::iterator start,
            vector<PeelOperation>::iterator end
        );

    vector<PeelOperation> all_possible_peels();

    bool build_peel_order();
    void print();
    bool peel(double *likelihood);
};

#endif

