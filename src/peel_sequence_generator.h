#ifndef LKG_PEELSEQUENCEGENERATOR_H_
#define LKG_PEELSEQUENCEGENERATOR_H_

using namespace std;

#include <vector>

#include "peeling.h"


class Pedigree;

// XXX i think I would need one that could peel :
// i) over genotypes at loci (to create L-sampler)
// ii) over descent graphs in-between loci (to calculate the LOD score)
//
// notes: peeling sequences might not be the same, because the number of
// things to peel over is different, ie: 2 alleles make 3 genotypes, but
// they make 4 possible descent state settings
// (if my understanding of this is correct, genotypes are unphased, but
// the descent graphs require phase)

class PeelSequenceGenerator {

    Pedigree& ped;
    vector<PeelOperation> peelorder;
    vector<PeelOperation> tmp;
    PeelingState state;

  public :
    PeelSequenceGenerator(Pedigree& p) 
        : ped(p), state(p) {}

    PeelOperation get_random_operation(vector<PeelOperation>& v);
    PeelOperation get_best_operation_heuristic(vector<PeelOperation>& v);
    PeelOperation get_best_operation(vector<PeelOperation>& v);
    
    void all_possible_peels(int& unpeeled);

    void build_peel_order();
    vector<PeelOperation>& get_peel_order();
    void print();
};

#endif

