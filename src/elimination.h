#ifndef LKG_GENOTYPEELIMINATION_H_
#define LKG_GENOTYPEELIMINATION_H_

using namespace std;

#include "pedigree.h"
#include "person.h"
#include "genotype.h"
#include "descent_graph.h"


class GenotypeElimination {
    
    Pedigree* ped;
    int **possible_genotypes;
    bool init_processing;
    bool sex_linked;
    
    void _init();
    void _copy(const GenotypeElimination& rhs);
    void _kill();
    void _initial_elimination();
    int _child_homoz(int** ds, unsigned locus, int mother, int father, 
                        int child, enum phased_genotype homoz, bool ismale);
    int _parent_homoz(int** ds, unsigned locus, int parent, int other_parent, 
                        int child, enum phased_genotype homoz, enum parentage p);
    bool _legal(int** ds, unsigned locus);
    bool _complete(int** ds, unsigned locus);
    bool _elimination_pass(int** ds, unsigned locus);
    enum parentage _state(enum phased_genotype child, 
                          enum phased_genotype parent, 
                          enum parentage p);
    void _random_eliminate(int** ds, unsigned locus);
    void _copy_descentstate(int** src, int** dst);
    void _copy_descentstate_locus(int** src, int** dst, unsigned locus);
    void _write_descentgraph(DescentGraph& d, int** ds);
    
 public :
    GenotypeElimination(Pedigree* p, bool sex_linked) :
        ped(p), 
        possible_genotypes(NULL),
        init_processing(false),
        sex_linked(sex_linked) {
        
        _init();
    }
    
    GenotypeElimination(const GenotypeElimination& rhs) :
        ped(rhs.ped),
        possible_genotypes(NULL),
        init_processing(rhs.init_processing),
        sex_linked(rhs.sex_linked) {
        
        _init();
        _copy(rhs);
    }
    
    virtual ~GenotypeElimination() {
        _kill();      
    }
    
    GenotypeElimination& operator=(const GenotypeElimination& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            init_processing = rhs.init_processing;
            sex_linked = rhs.sex_linked;
        
            _kill();
            _init();
            _copy(rhs);
        }
        
        return *this;
    }
    
    bool elimination();
    bool random_descentgraph(DescentGraph& d);
    bool is_legal(int id, int locus, int value);
};

#endif

