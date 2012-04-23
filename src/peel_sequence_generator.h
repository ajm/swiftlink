#ifndef LKG_PEELSEQUENCEGENERATOR_H_
#define LKG_PEELSEQUENCEGENERATOR_H_

using namespace std;

#include <vector>

#include "peeling.h"
#include "elimination.h"


class Pedigree;
class GeneticMap;

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

    Pedigree* ped;
    GeneticMap* map;
    bool verbose;
    vector<PeelOperation> peelorder;
    vector<PeelOperation> tmp;
    PeelingState state;
    GenotypeElimination ge;
    
    void find_previous_functions(PeelOperation& op);
    void find_generic_functions(PeelOperation& op);
    void find_child_functions(PeelOperation& op);
    int  find_function_containing(vector<unsigned>& nodes);
    void bruteforce_assignments(PeelOperation& op);
    PeelOperation get_random_operation(vector<PeelOperation>& v);
    PeelOperation get_best_operation_heuristic(vector<PeelOperation>& v);
    PeelOperation get_best_operation(vector<PeelOperation>& v);
    void all_possible_peels(int& unpeeled);

    int calculate_cost(vector<unsigned int>& seq);
    void rebuild_peel_order(vector<unsigned int>& seq);

  public :
    PeelSequenceGenerator(Pedigree* p, GeneticMap* m, bool verbose) : 
        ped(p),
        map(m),
        verbose(verbose),
        peelorder(),
        tmp(),
        state(p),
        ge(p) {
        
        ge.elimination();    
    }
        
    ~PeelSequenceGenerator() {}
    
    PeelSequenceGenerator(const PeelSequenceGenerator& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        verbose(rhs.verbose),
        peelorder(rhs.peelorder),
        tmp(rhs.tmp),
        state(rhs.state),
        ge(rhs.ge) {}
        
    PeelSequenceGenerator& operator=(const PeelSequenceGenerator& rhs) {
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            verbose = rhs.verbose;
            peelorder = rhs.peelorder;
            tmp = rhs.tmp;
            state = rhs.state;
            ge = rhs.ge;
        }
        
        return *this;
    }
    
    void build_peel_order();
    bool read_from_file(string filename);
    vector<PeelOperation>& get_peel_order();
    
    unsigned score_peel_sequence();
    
    string debug_string();
};

#endif
