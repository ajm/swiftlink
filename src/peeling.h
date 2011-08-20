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
    unsigned int peelnode;       // what is being peeled by this operation
    bool used;
    
 public :
    PeelOperation() :  
        type(NULL_PEEL), 
        cutset(), 
        peelnode(0), 
        used(false) {}
        
    ~PeelOperation() {}
    
    void set_used() {
        used = true;
    }
    
    bool is_used() const {
        return used;
    }
    
    bool in_cutset(unsigned node) const {
        for(unsigned i = 0; i < cutset.size(); ++i) {
            if(cutset[i] == node) {
                return true;
            }
        }
        return false;
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
    
    void add_cutnode(unsigned int c) {
        for(unsigned int i = 0; i < cutset.size(); ++i) {
            if(cutset[i] == c)
                return;
        }
        cutset.push_back(c);
    }
    
    void remove_cutnode(unsigned int c) {
        vector<unsigned int>::iterator it = find(cutset.begin(), cutset.end(), c);
        if(it != cutset.end())
            cutset.erase(it);
    }
    
    unsigned get_cutnode(unsigned i) const {
        return cutset[i];
    }
    
    unsigned get_peelnode() const {
        return peelnode;
    }
    
    void set_peelnode(unsigned i) {
        peelnode = i;
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
        
        ss  << peeloperation_str(type) << " " \
            << "peelnode = " << peelnode << " " \
            << "cutset = (";
        
        unsigned tmp = cutset.size();
        for(unsigned i = 0; i < tmp; ++i) {
            ss << cutset[i];
            if(i != (tmp-1)) {
                ss << ",";
            }
        }
        ss << ") ";
        
        return ss.str();
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

