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


class SimpleGraph {
    unsigned int id;
    vector<unsigned int> cutset;
    
 public:
    SimpleGraph(unsigned int i) :
        id(i),
        cutset() {}
        
    ~SimpleGraph() {}
        
    SimpleGraph(const SimpleGraph& rhs) :
        id(rhs.id),
        cutset(rhs.cutset) {}
        
    SimpleGraph& operator=(const SimpleGraph& rhs) {
        if(&rhs != this) {
            id = rhs.id;
            cutset = rhs.cutset;
        }
        
        return *this;
    }
    
    bool operator<(const SimpleGraph& rhs) {
        return get_cutset_size() < rhs.get_cutset_size();
    }
    
    unsigned int get_id() const { 
        return id;
    }
    
    void add(unsigned int i) {
        //if(not binary_search(cutset.begin(), cutset.end(), i)) {
        if(find(cutset.begin(), cutset.end(), i) == cutset.end()) {
            cutset.push_back(i);
            //sort(cutset.begin(), cutset.end());
        }
    }
    
    void remove_node(unsigned int node) {        
        vector<unsigned int>::iterator it;
        
        it = find(cutset.begin(), cutset.end(), node);
        
        if(it != cutset.end())
            cutset.erase(it);
    }
    
    void reset() {
        cutset.clear();
    }
    
    int get_cutset_size() const {
        return cutset.size();
    }
    
    int get_cost() const {
        return static_cast<int>(pow(4.0, static_cast<double>(cutset.size())));
    }
    
    vector<unsigned int>& get_cutset() {
        return cutset;
    }
    
    string debug_string() {
        stringstream ss;
        unsigned tmp;
        
        ss << "peelnode = " << id << " " \
           << "cutset = (";
        
        tmp = get_cutset_size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << cutset[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        return ss.str();
    }
};

class PeelSequenceGenerator {

    Pedigree* ped;
    GeneticMap* map;
    bool verbose;
    vector<PeelOperation> peelorder;
    PeelingState state;
    GenotypeElimination ge;
    
    vector<SimpleGraph> graph;
    
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
    
    
    void build_simple_graph();
    void build_peel_sequence();
    void eliminate_node(vector<SimpleGraph>& tmp, unsigned int node);
    unsigned int get_cost(vector<unsigned int>& peel);
    unsigned int get_proper_cost(vector<unsigned int>& peel);
    void print_graph(vector<SimpleGraph>& g);
    void find_prev_functions(PeelOperation& op);

  public :
    PeelSequenceGenerator(Pedigree* p, GeneticMap* m, bool verbose) : 
        ped(p),
        map(m),
        verbose(verbose),
        peelorder(),
        state(p),
        ge(p),
        graph() {
        
        ge.elimination();
        build_simple_graph();
        //print_graph(graph);
        build_peel_sequence();
    }
        
    ~PeelSequenceGenerator() {}
    
    PeelSequenceGenerator(const PeelSequenceGenerator& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        verbose(rhs.verbose),
        peelorder(rhs.peelorder),
        state(rhs.state),
        ge(rhs.ge),
        graph(rhs.graph) {}
        
    PeelSequenceGenerator& operator=(const PeelSequenceGenerator& rhs) {
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            verbose = rhs.verbose;
            peelorder = rhs.peelorder;
            state = rhs.state;
            ge = rhs.ge;
            graph = rhs.graph;
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
