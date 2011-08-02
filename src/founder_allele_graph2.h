#ifndef LKG_FOUNDERALLELEGRAPH_H_
#define LKG_FOUNDERALLELEGRAPH_H_

using namespace std;

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include "genotype.h"
#include "pedigree.h"


class DescentGraph;
class GeneticMap;

class FounderAlleleNode {
    int id;
    enum unphased_genotype label;
    
    string debug_string() {
        stringstream ss;
        ss << "{" << id << ", " << genotype_string(label) << "}";
        return ss.str();
    }
};

class AdjacencyRecord {
    vector<FounderAlleleNode> edges;
    
 public:
    AdjacencyRecord() : 
        edges() {}
        
    AdjacencyRecord(const AdjacencyRecord& rhs) :
        edges(rhs.edges) {}
    
    AdjacencyRecord& operator=(const AdjacencyRecord& rhs) {
        if(this != &rhs) {
            edges = rhs.edges;
        }
        return *this;
    }
    
    bool add(FounderAlleleNode& f) {
        for(unsigned i = 0; i < edges.size(); ++i) {
            if(edges[i].id == f.id) {
                return edges[i].label == f.label;
            }
        }
        
        edges.push_back(f);
        
        return true;
    }
    
    void remove() {}
    
    unsigned size() {
        return edges.size();
    }
    
    string debug_string() {
        stringstream ss;
        for(unsigned i = 0; i < edges.size(); ++i) {
            ss << edges[i].debug_string() << " ";
        }
        return ss.str();
    }
};

class AdjacencyMatrix {
    vector<AdjacencyRecord> nodes;
    
 public:
    AdjacencyMatrix(int num_nodes) :
        nodes(num_nodes) {}
        
    AdjacencyMatrix(const AdjacencyMatrix& rhs) :
        nodes(rhs.nodes) {}
        
    AdjacencyMatrix& operator=(const AdjacencyMatrix& rhs) {
        if(this != &rhs) {
            nodes = rhs.nodes;
        }
        return *this;
    }
    
    const AdjacencyRecord& operator[](unsigned i) const {
        return nodes[i];
    }
        
    bool add(int founderallele1, int founderallele2, enum unphased_genotype g) { 
        if(g == UNTYPED)
            return true;
    
        if(not nodes[founderallele1].add(FounderAlleleNode(founderallele2, g)))
            return false;
        
        if(founderallele1 != founderallele2) {
            if(not nodes[founderallele2].add(FounderAlleleNode(founderallele1, g)))
                return false;
        }
        
        return true;
    }
    
    void remove() {}
    
    unsigned size() {
        return nodes.size();
    }
    
    string debug_string() {
        stringstream ss;
        for(unsigned i = 0; i < nodes.size(); ++i) {
            ss << i << ": " << nodes[i].debug_string() << "\n";
        }
        return ss.str();
    }
};

class GraphComponent {
    vector<int> founderalleles;
    double prob;
    
 public:
    GraphComponent() : 
        founderalleles(), 
        prob(0.0) {}
    
    ~GraphComponent() {}
    
    void add(int i) {
        if(not contains(i)) {
            founderalleles.push_back(i);
        }
    }
    
    void remove(int i) {
        vector<int>::iterator it = find(founderalleles.begin(), founderalleles.end(), i);
        
        if(it != founderalleles.end()) {
            founderalleles.erase(it);
        }
    }
    
    bool contains(int i) {
        return find(founderalleles.begin(), founderalleles.end(), i) != founderalleles.end();
    } 
};

class FounderAlleleGraph {
	
    Pedigree* ped;
    GeneticMap* map;
    int num_alleles;
    
    AdjacencyMatrix matrix;
    vector<GraphComponent> components;
    
    
    bool _check_legality(int node, int node_assignment, vector<int>& assignments);
    void _assign_and_recurse(int *component, int component_size, unsigned locus, int current_index, vector<int>& assignment, double *prob, double *best);
    double _enumerate_component(int *component, int component_size, unsigned locus);
    
 public :
	FounderAlleleGraph(Pedigree* ped, GeneticMap* map) :
	    ped(ped), 
	    map(map),
	    num_alleles(2 * p->num_founders()),
	    matrix(num_founder_alleles),
	    components() {}
	    
	FounderAlleleGraph(const FounderAlleleGraph& rhs) :
	    ped(rhs.ped),
	    map(rhs.map),
	    num_alleles(rhs.num_alleles),
	    matrix(rhs.matrix),
	    components(rhs.components) {}
	
	~FounderAlleleGraph() {}
	
	FounderAlleleGraph& operator=(const FounderAlleleGraph& rhs) {
	    if(this != &rhs) {
	        ped = rhs.ped;
	        map = rhs.map;
	        num_alleles = rhs.num_alleles;
	        matrix = rhs.matrix;
	        components = rhs.components;
	    }
	    return *this;
	}
    
    string debug_string() {
        return matrix.debug_string();
    }
    
    bool populate(DescentGraph& d, unsigned locus);
	bool likelihood(double* probs, unsigned locus);
};

#endif

