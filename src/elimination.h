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
    
    void _initial_elimination();
    int _child_homoz(int** ds, unsigned locus, int mother, int father, 
                        int child, enum phased_genotype homoz);
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
    GenotypeElimination(Pedigree* p) : ped(p), init_processing(false) {
        possible_genotypes = new int*[ped->num_markers()];
        
        for(int i = 0; i < int(ped->num_markers()); ++i) {
            possible_genotypes[i] = new int[ped->num_members()];
        }
    }
    
    virtual ~GenotypeElimination() {
        for(int i = 0; i < int(ped->num_markers()); ++i) {
            delete[] possible_genotypes[i];
        }
        
        delete[] possible_genotypes;        
    }
    
    bool elimination();
    bool random_descentgraph(DescentGraph& d);
//    virtual void random_descentgraph(DescentGraph* d) = 0;
};

#endif
