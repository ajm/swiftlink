#ifndef LKG_PEELING_H_
#define LKG_PEELING_H_

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include "pedigree.h"


enum peeloperation {
    NULL_PEEL,
    CHILD_PEEL,
    PARENT_PEEL,
    PARTNER_PEEL,
    LAST_PEEL
};

class PeelOperation {
    enum peeloperation type;
    vector<unsigned int> cutset; // what is being peeled on to by this operation
    vector<unsigned int> nuclearset; // what is the nuclear family of the peelnode that is contained in the cutset
    unsigned int peelnode;       // what is being peeled by this operation
    bool used;
    int prev1;
    int prev2;
    
    vector<vector<int> > assignments;
    
 public :
    PeelOperation(unsigned int peelnode) :  
        type(NULL_PEEL), 
        cutset(), 
        nuclearset(),
        peelnode(peelnode), 
        used(false),
        prev1(-1),
        prev2(-1),
        assignments() {
    
        nuclearset.push_back(peelnode);    
    }
    
    PeelOperation(const PeelOperation& rhs) :
        type(rhs.type),
        cutset(rhs.cutset),
        nuclearset(rhs.nuclearset),
        peelnode(rhs.peelnode),
        used(rhs.used),
        prev1(rhs.prev1),
        prev2(rhs.prev2),
        assignments(rhs.assignments) {}
        
    ~PeelOperation() {}
    
    void set_used() {
        used = true;
    }
    
    bool is_used() const {
        return used;
    }
    
    bool in_cutset(unsigned node) const {
        return find(cutset.begin(), cutset.end(), node) != cutset.end();
    }
    
    void set_type(enum peeloperation po) {
        type = po;
    }
    
    enum peeloperation get_type() const {
        return type;
    }
    
    unsigned int get_cutset_size() const { 
        return cutset.size();
    }
    
    vector<unsigned int>& get_cutset() {
        return cutset;
    }
    
    unsigned int get_nuclearset_size() {
        return nuclearset.size();
    }
    
    vector<unsigned int>& get_nuclearset() {
        return nuclearset;
    }
    
    void add_cutnode(unsigned int c, bool in_nuclear_fam) {
        if(not in_cutset(c)) {
            cutset.push_back(c);
            
            if(in_nuclear_fam) {
                nuclearset.push_back(c);
            }
        }
    }
    
    void remove_cutnode(unsigned int c) {
        vector<unsigned int>::iterator it = find(cutset.begin(), cutset.end(), c);
        if(it != cutset.end())
            cutset.erase(it);
    }
    
    unsigned get_cutnode(unsigned int i) const {
        return cutset[i];
    }
    
    bool contains_cutnodes(vector<unsigned>& nodes) {
        for(unsigned i = 0; i < nodes.size(); ++i) {
            if(find(cutset.begin(), cutset.end(), nodes[i]) == cutset.end())
                return false;
        }
        
        return true;
    }
    
    unsigned get_peelnode() const {
        return peelnode;
    }
    
    void set_peelnode(unsigned i) {
        peelnode = i;
    }
    
    void set_previous_operation(int i) {
        if(prev1 == -1) {
            prev1 = i;
            return;
        }
        
        if(prev2 == -1) {
            prev2 = i;
            return;
        }
        
        abort();
    }
    
    int get_previous_op1() const {
        return prev1;
    }

    int get_previous_op2() const {
        return prev2;
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
        
        /*
        ss << "nuclear = (";
        tmp = nuclearset.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << nuclearset[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        */
        
        ss << " prev = (" << prev1 << ", " << prev2 << ")";
        
        return ss.str();
    }
    
    void set_index_values(vector<vector<int> > assigns) {
        assignments = assigns;
    }
    
    vector<vector<int> >* get_index_values() {
        return &assignments;
    }
    
    bool operator<(const PeelOperation& p) const {
		return get_cutset_size() < p.get_cutset_size();
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
    
    string debug_string() {
        stringstream ss;
        
        for(unsigned i = 0; i < peeled.size(); ++i) {
            ss << i << "\t" << (peeled[i] ? "peeled" : "unpeeled") << "\n";
        }
        
        return ss.str();
    }
};

#endif

