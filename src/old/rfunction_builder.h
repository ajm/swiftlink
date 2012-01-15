#ifndef LKG_RFUNCTIONBUILDER_H_
#define LKG_RFUNCTIONBUILDER_H_

#include <vector>

#include "peeling.h"
#include "person.h"
#include "pedigree.h"
#include "rfunction.h"


class Pedigree;
class GeneticMap;

template<class T>
class RfunctionBuilder {
    
    Pedigree* ped;
    GeneticMap* map;
    vector<T*>& functions;

/*
    void find_previous_functions(const PeelOperation& peel, Rfunction** f1, Rfunction** f2) {
        switch(peel.get_type()) {
            case CHILD_PEEL:
                find_child_functions(peel, f1, f2);
                break;
                
            case PARTNER_PEEL:
            case PARENT_PEEL:
            case LAST_PEEL:
                find_generic_functions(peel, f1, f2);
                break;
                
            default:
                fprintf(stderr, "error: default should never be reached! (%s:%d)\n", __FILE__, __LINE__);
                abort();
        }
    }
    
    void find_generic_functions(const PeelOperation& peel, Rfunction** f1, Rfunction** f2) {
        vector<unsigned> tmp;
        tmp.push_back(peel.get_peelnode());
        
        *f1 = find_function_containing(tmp);
        *f2 = find_function_containing(tmp);
    }
    
    void find_child_functions(const PeelOperation& peel, Rfunction** f1, Rfunction** f2) {
        Person* p = ped->get_by_index(peel.get_peelnode());
        vector<unsigned> tmp;
        tmp.push_back(p->get_maternalid());
        tmp.push_back(p->get_paternalid());
        
        *f1 = find_function_containing(tmp);
        
        // don't even bother looking if the child is a leaf
        if(p->isleaf()) {
            return;
        }
        
        tmp.clear();
        tmp.push_back(peel.get_peelnode());
        
        *f2 = find_function_containing(tmp);
    }
    
    Rfunction* find_function_containing(vector<unsigned>& nodes) {
        T* ret = NULL;
        
        for(unsigned i = 0; i < functions.size(); ++i) {
            if(functions[i]->is_used()) {
                continue;
            }
            
            if(functions[i]->contains_cutnodes(nodes)) {
                ret = functions[i];
                functions[i]->set_used();
                break;
            }
        }
        
        return static_cast<Rfunction*>(ret);
    }
*/

 public :
    RfunctionBuilder(Pedigree* p, GeneticMap* m, vector<T*>& func) : 
        ped(p), map(m), functions(func) {}
    
    T* createRfunction(const PeelOperation& po) {
        Rfunction* prev1 = po.get_previous_op1() == -1 ? NULL : functions[po.get_previous_op1()];
        Rfunction* prev2 = po.get_previous_op2() == -1 ? NULL : functions[po.get_previous_op2()];
        
        //find_previous_functions(po, &prev1, &prev2);
        
        return new T(po, ped, map, prev1, prev2);
    }
};

#endif
