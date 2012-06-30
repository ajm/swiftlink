#ifndef LKG_PEELING_H_
#define LKG_PEELING_H_

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>

#include "pedigree.h"


enum peeloperation {
    NULL_PEEL,
    CHILD_PEEL,
    PARENT_PEEL,
    PARTNER_PEEL,
    LAST_PEEL
};

class PeelOperation {
    
    Pedigree* ped;
    enum peeloperation type;        // which function is called later
    unsigned int peelnode;          // what is being peeled by this operation
    bool used;                      // book-keeping
    vector<unsigned int> cutset;    // what is being peeled on to by this operation
    vector<unsigned int> children;  // used to know what to pre-calculate for transmission probs (parent + child peels)
    vector<unsigned int> previous;  // all previous functions (0-3) (at least 3 is the most i have seen...)
    
    // cached indices 
    vector<vector<int> > assignments;
    vector<vector<int> > presum_indices;
    vector<vector<int> > matrix_indices;
    vector<int> lod_indices;
    
    
 public :
    PeelOperation(Pedigree* ped, unsigned int peelnode) :
        ped(ped),
        type(NULL_PEEL), 
        peelnode(peelnode),
        used(false), 
        cutset(), 
        children(),
        previous(),
        assignments(),
        presum_indices(),
        matrix_indices(),
        lod_indices() {}
    
    PeelOperation(const PeelOperation& rhs) :
        ped(rhs.ped),
        type(rhs.type),
        peelnode(rhs.peelnode),
        used(rhs.used),
        cutset(rhs.cutset),
        children(rhs.children),
        previous(rhs.previous),
        assignments(rhs.assignments),
        presum_indices(rhs.presum_indices),
        matrix_indices(rhs.matrix_indices),
        lod_indices(rhs.lod_indices) {}
        
    ~PeelOperation() {}
    
    PeelOperation& operator=(const PeelOperation& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            type = rhs.type;
            peelnode = rhs.peelnode;
            used = rhs.used;
            cutset = rhs.cutset;
            children = rhs.children;
            previous = rhs.previous;
            assignments = rhs.assignments;
            presum_indices = rhs.presum_indices;
            matrix_indices = rhs.matrix_indices;
            lod_indices = rhs.lod_indices;
        }
        
        return *this;
    }
    
    bool operator<(const PeelOperation& rhs) const {
		return get_cutset_size() < rhs.get_cutset_size();
	}
    
    unsigned int get_peelnode() const {
        return peelnode;
    }
    
    void set_used() {
        used = true;
    }
    
    bool is_used() const {
        return used;
    }
    
    bool in_cutset(unsigned int node) const {
        return find(cutset.begin(), cutset.end(), node) != cutset.end();
    }
    
    void add_cutnode(unsigned int c) {
        if(not in_cutset(c)) {
            cutset.push_back(c);
        }
    }
    
    void remove_cutnode(unsigned int c) {
        vector<unsigned int>::iterator it = find(cutset.begin(), cutset.end(), c);
        
        if(it != cutset.end()) {
            cutset.erase(it);
        }
    }
    
    unsigned get_cutnode(unsigned int i) const {
        return cutset[i];
    }
    
    void reset() {
        cutset.clear();
    }
    
    unsigned int get_cost() const {
        return static_cast<int>(pow(4.0, static_cast<double>(cutset.size())));
    }
    
    bool contains_cutnodes(vector<unsigned>& nodes) {
        for(unsigned int i = 0; i < cutset.size(); ++i) {
            if(find(nodes.begin(), nodes.end(), cutset[i]) == nodes.end())
                return false;
        }
        
        return true;
    }
    
    unsigned int get_cutset_size() const { 
        return cutset.size();
    }
    
    vector<unsigned int>& get_cutset() {
        return cutset;
    }
    
    void set_type(enum peeloperation po) {
        type = po;
        
        if(type == PARENT_PEEL) {
            Person* p = ped->get_by_index(peelnode);
            
            for(unsigned int i = 0; i < cutset.size(); ++i) {
                if(p->is_offspring(cutset[i])) {
                    children.push_back(cutset[i]);
                }
            }
            
            /*
            if(not (state.is_peeled(p->get_maternalid()) or state.is_peeled(p->get_paternalid())) {
                children.push_back(peelnode);
            }
            */
        }
    }
    
    enum peeloperation get_type() const {
        return type;
    }
    
    unsigned int get_children_size() const {
        return children.size();
    }
    
    vector<unsigned int>& get_children() {
        return children;
    }
    
    void add_prevfunction(unsigned int i) {
        previous.push_back(i);
    }
    
    vector<unsigned int>& get_prevfunctions() {
        return previous;
    }
    
    unsigned int get_prev_size() const {
        return previous.size();
    }

    string peeloperation_str(enum peeloperation po) {
        
        switch(po) {
            case NULL_PEEL:
                return "null";
            case CHILD_PEEL:
                return "child";
            case PARENT_PEEL:
                return "parent";
            case PARTNER_PEEL:
                return "partner";
            case LAST_PEEL:
                return "last";
        }
        
        abort();
    }
    
    string debug_string() {
        stringstream ss;
        unsigned tmp;
        
        ss << peeloperation_str(type) << " " \
           << "peelnode = " << peelnode << " " \
           << "cutset = (";
        
        tmp = cutset.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << cutset[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        ss << "prev = (";
        
        tmp = previous.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << previous[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        ss << " children = (";
        tmp = children.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << children[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        return ss.str();
    }
    
    string translated_debug_string(Pedigree* ped) {
        stringstream ss;
        unsigned tmp;
        
        ss << peeloperation_str(type) << "\t" \
           << "peelnode = " << peelnode << "\t" \
           << "id = " << ped->get_by_index(peelnode)->get_id() << "\t" \
           << "cutset = (";
        
        tmp = cutset.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << ped->get_by_index(cutset[i])->get_id();
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ")\t";
        
        ss << "prev = (";
        
        tmp = previous.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << previous[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        ss << " children = (";
        tmp = children.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << ped->get_by_index(children[i])->get_id();
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        return ss.str();
    }
    
    string code_output() {
        stringstream ss;
        
        ss << "current.push_back(" << peelnode << ");";
        
        return ss.str();
    }
    
    void set_index_values(vector<vector<int> > assigns) {
        assignments = assigns;
    }
    
    void set_presum_indices(vector<vector<int> > indices) {
        presum_indices = indices;
    }
    
    void set_matrix_indices(vector<vector<int> > indices) {
        matrix_indices = indices;
    }
    
    void set_lod_indices(vector<int> indices) {
        lod_indices = indices;
    }
    
    // return-by-value as each locus uses it to store temporary values
    // for the node being peeled as well 
    vector<vector<int> > get_index_values() {
        return assignments;
    }
    
    vector<int>* get_presum_indices(int locus) {
        return &presum_indices[locus];
    }
    
    vector<int>* get_matrix_indices(int locus) {
        return &matrix_indices[locus];
    }
    
    vector<int>* get_lod_indices() {
        return &lod_indices;
    }
};

class PeelingState {
    vector<bool> peeled;

  public :
    PeelingState(Pedigree* p) : 
        peeled(p->num_members(), false) {}

    bool is_peeled(unsigned int i) {
        return peeled[i];
    }

    void set_peeled(unsigned int i) {
        peeled[i] = true;
    }
    
    void toggle_peeled(unsigned int i) {
        peeled[i] = peeled[i] ? false : true;
    }
    
    void toggle_peel_operation(PeelOperation& operation) {
        toggle_peeled(operation.get_peelnode());
    }
    
    bool is_final_node(unsigned int node) {
        for(unsigned int i = 0; i < peeled.size(); ++i) {
            if((i != node) and not peeled[i])
                return false;
        }
        
        return true;
    }
    
    void reset() {
        for(unsigned i = 0; i < peeled.size(); ++i)
            peeled[i] = false;
    }
    
    string debug_string() {
        stringstream ss;
        
        for(unsigned i = 0; i < peeled.size(); ++i) {
            ss << i << "\t" << (peeled[i] ? "peeled" : "unpeeled") << "\n";
        }
        
        return ss.str();
    }
};

#endif

