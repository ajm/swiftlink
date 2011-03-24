#ifndef LKG_PEELING_H_
#define LKG_PEELING_H_

using namespace std;

#include <vector>
#include <cstdio>

#include "pedigree.h"


enum peeloperation {
    NULL_PEEL,
    CHILD_PEEL,
    PARENT_PEEL,
    PARTNER_PEEL,
    LAST_PEEL
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
};

class PeelOperation {
    int pivot;
    enum peeloperation type;
    vector<unsigned int> cutset;
    
  public :
    PeelOperation() 
        : pivot(-1) {}
    
    PeelOperation(unsigned int pivot_node) 
        : pivot(pivot_node) {}

    ~PeelOperation() {}
    
    unsigned int get_pivot() const { 
        return pivot;
    }

    void set_pivot(unsigned int p) {
        pivot = p;
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

    void set_type(enum peeloperation po) {
        type = po;
    }

    enum peeloperation get_type() {
        return type;
    }
    
    void print() const {
        unsigned int tmp = cutset.size();

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

        printf("pivot = %d, cutset = (", pivot);
        for(unsigned int i = 0; i < tmp; ++i) {
            printf("%d", cutset[i]);
            if(i != (tmp-1)) {
                putchar(',');
            }
        }
        printf(")\n");
    }

    bool operator<(const PeelOperation& p) const {
		return cutset.size() < p.cutset.size();
	}
};

#endif

