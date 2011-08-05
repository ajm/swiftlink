#ifndef LKG_FounderAlleleGraph2_H_
#define LKG_FounderAlleleGraph2_H_

using namespace std;

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include "genotype.h"
#include "pedigree.h"


class DescentGraph;
class GeneticMap;

class FounderAlleleNode;
class AdjacencyRecord;
class AdjacencyMatrix;
class GraphComponent;

class FounderAlleleNode {

 public:
    int id;
    enum unphased_genotype label;
    
    FounderAlleleNode(int id, enum unphased_genotype label) :
        id(id),
        label(label) {}
    
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
    
    const FounderAlleleNode& operator[](unsigned i) const {
        return edges[i];
    }
    
    bool add(FounderAlleleNode& f) {
        //fprintf(stderr, "edges.size() = %d\n", int(edges.size()));
        
        for(unsigned i = 0; i < edges.size(); ++i) {
            if(edges[i].id == f.id) {
                return edges[i].label == f.label;
            }
        }
        
        edges.push_back(f);
        
        return true;
    }
    
    void remove() {}
    
    unsigned size() const {
        return edges.size();
    }
    
    void clear() {
        edges.clear();
    }
    
    bool contains(int i, enum unphased_genotype* g) {
        for(unsigned j = 0; j < edges.size(); ++j) {
            if(edges[j].id == i) {
                *g = edges[j].label;
                return true;
            }
        }
        return false;
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
    
        FounderAlleleNode f1(founderallele2, g);
    
        if(not nodes[founderallele1].add(f1))
            return false;
        
        if(founderallele1 != founderallele2) {
            FounderAlleleNode f2(founderallele1, g);
            
            if(not nodes[founderallele2].add(f2))
                return false;
        }
        
        return true;
    }
    
    void remove() {}
    
    unsigned size() const {
        return nodes.size();
    }
    
    void clear() {
        // we only want to clear the internal structure,
        // all the rest of the code relies on each founder
        // allele having at least an empty reference in
        // the nodes vector
        for(unsigned i = 0; i < nodes.size(); ++i) {
            nodes[i].clear();
        }
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
    
    const int& operator[](unsigned i) {
        return founderalleles[i];
    } 
    
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
    
    void clear() {
        founderalleles.clear();
    }
    
    unsigned size() {
        return founderalleles.size();
    }
};

class FounderAlleleGraph2 {
	
    Pedigree* ped;
    GeneticMap* map;
    unsigned num_alleles;
    unsigned locus;
    
    AdjacencyMatrix matrix;
    vector<GraphComponent> components;
    
    
    bool populate_graph(DescentGraph& d);
    void populate_components();
    bool correct_alleles(enum unphased_genotype g, int allele1);
    bool correct_alleles(enum unphased_genotype g, int allele1, int allele2);
    bool legal(GraphComponent& gc, vector<unsigned>& assignment);
    double component_likelihood(vector<unsigned>& q);
    double enumerate_component(GraphComponent& c);
    
 public :
	FounderAlleleGraph2(Pedigree* ped, GeneticMap* map, unsigned locus) :
	    ped(ped), 
	    map(map),
	    num_alleles(2 * ped->num_founders()),
	    locus(locus),
	    matrix(num_alleles),
	    components() {}
	    
	FounderAlleleGraph2(const FounderAlleleGraph2& rhs) :
	    ped(rhs.ped),
	    map(rhs.map),
	    num_alleles(rhs.num_alleles),
	    locus(rhs.locus),
	    matrix(rhs.matrix),
	    components(rhs.components) {}
	
	~FounderAlleleGraph2() {}
	
	FounderAlleleGraph2& operator=(const FounderAlleleGraph2& rhs) {
	    if(this != &rhs) {
	        ped = rhs.ped;
	        map = rhs.map;
	        num_alleles = rhs.num_alleles;
	        locus = rhs.locus;
	        matrix = rhs.matrix;
	        components = rhs.components;
	    }
	    return *this;
	}
    
    string debug_string() {
        return matrix.debug_string();
    }
    
    // XXX this is temporary, I want to have a graph per locus in the future
    void set_locus(unsigned l) {
        locus = l;
    }
    
    void reset() {
        matrix.clear();
        components.clear();
    }
    
    bool populate(DescentGraph& d);
	bool likelihood(double* probs);
};

#endif

