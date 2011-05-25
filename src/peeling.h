#ifndef LKG_PEELING_H_
#define LKG_PEELING_H_

using namespace std;

#include <vector>
#include <algorithm>
#include <cstdio>

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
    vector<unsigned int> peelset; // what is being peeled by this operation
    
    public :
    PeelOperation() {}
    ~PeelOperation() {}
    
    void set_type(enum peeloperation po) {
        type = po;
    }
    
    enum peeloperation get_type() {
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
    
    unsigned int get_peelset_size() const { 
        return peelset.size();
    }
    
    vector<unsigned int>& get_peelset() {
        return peelset;
    }
    
    unsigned get_peelnode(unsigned i) {
        return peelset[i];
    }
    
    unsigned get_cutnode(unsigned i) {
        return cutset[i];
    }
    
    void add_peelnode(unsigned int c) {
        for(unsigned int i = 0; i < peelset.size(); ++i) {
            if(peelset[i] == c)
                return;
        }
        peelset.push_back(c);
    }
    
    void print() const {
        unsigned int tmp;
        
        switch(type) {
            case NULL_PEEL:
                printf("null ");
                break;
            case CHILD_PEEL:
                printf("child ");
                break;
            case PARENT_PEEL:
                printf("parent ");
                break;
            case PARTNER_PEEL:
                printf("partner ");
                break;
            case LAST_PEEL:
                printf("last ");
                break;
            default:
                printf("error ");
                break;
        }
        
        printf("peelset = (");
        tmp = peelset.size();
        for(unsigned i = 0; i < tmp; ++i) {
            printf("%d", peelset[i]);
            if(i != (tmp-1)) {
                putchar(',');
            }
        }
        printf(") ");
        
        printf("cutset = (");
        tmp = cutset.size();
        for(unsigned i = 0; i < tmp; ++i) {
            printf("%d", cutset[i]);
            if(i != (tmp-1)) {
                putchar(',');
            }
        }
        printf(") ");
        
        printf("\n");
        
    }
    
    bool operator<(const PeelOperation& p) const {
		return get_cutset_size() < p.get_cutset_size();
	}
};

class PeelingState {
    vector<bool> peeled;

  public :
    PeelingState(Pedigree& p) 
        : peeled(p.num_members(), false) {}

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
        vector<unsigned>& tmp = operation.get_peelset();
        
        for(unsigned i = 0; i < tmp.size(); ++i) {
            toggle_peeled(tmp[i]);
        }
    }
    
    void print() {
        for(unsigned i = 0; i < peeled.size(); ++i) {
            printf("%d\t%s\n", i, peeled[i] ? "peeled" : "unpeeled");
        }
    }
};

#endif

