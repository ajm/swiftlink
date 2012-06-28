#ifndef LKG_PEELSEQUENCEGENERATOR_H_
#define LKG_PEELSEQUENCEGENERATOR_H_

using namespace std;

#include <vector>
#include <algorithm>
#include <cmath>

#include "peeling.h"
#include "elimination.h"


class Pedigree;
class GeneticMap;


class PeelSequenceGenerator {

    Pedigree* ped;
    GeneticMap* map;
    bool verbose;
    vector<PeelOperation> peelorder;
    PeelingState state;
    GenotypeElimination ge;
    
    
    
    void find_prev_functions(PeelOperation& op);
    int find_function_containing(vector<unsigned>& nodes);
    void bruteforce_assignments(PeelOperation& op);
    
    int calculate_cost(vector<unsigned int>& seq);
    void finalise_peel_order(vector<unsigned int>& seq);
    
    void set_type(PeelOperation& p);
    
    
    void build_simple_graph();
    void build_peel_sequence();
    void eliminate_node(vector<PeelOperation>& tmp, unsigned int node);
    unsigned int get_cost(vector<unsigned int>& peel);
    unsigned int get_proper_cost(vector<unsigned int>& peel);
    bool is_legit(vector<unsigned int>& peel);
    
    void random_downhill_search(vector<unsigned int>& current);

  public :
    PeelSequenceGenerator(Pedigree* p, GeneticMap* m, bool verbose) : 
        ped(p),
        map(m),
        verbose(verbose),
        peelorder(),
        state(p),
        ge(p) {
        
        ge.elimination();
        
        build_simple_graph();
        build_peel_sequence();
    }
        
    ~PeelSequenceGenerator() {}
    
    PeelSequenceGenerator(const PeelSequenceGenerator& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        verbose(rhs.verbose),
        peelorder(rhs.peelorder),
        state(rhs.state),
        ge(rhs.ge) {}
        
    PeelSequenceGenerator& operator=(const PeelSequenceGenerator& rhs) {
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            verbose = rhs.verbose;
            peelorder = rhs.peelorder;
            state = rhs.state;
            ge = rhs.ge;
        }
        
        return *this;
    }
    
    //bool read_from_file(string filename);
    vector<PeelOperation>& get_peel_order();
    unsigned int get_peeling_cost();
    
    string debug_string();
};

#endif

