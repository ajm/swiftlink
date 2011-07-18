#ifndef LKG_TRAITRFUNCTION_H_
#define LKG_TRAITRFUNCTION_H_

#include <vector>

#include "trait.h"
#include "rfunction.h"
#include "descent_graph_types.h"

class DescentGraph;
class Pedigree;
class GeneticMap;

class TraitRfunction : public Rfunction {
    
    double get_recombination_probability(DescentGraph* dg, 
                                         unsigned locus,
                                         unsigned person_id, 
                                         int maternal_allele, 
                                         int paternal_allele);
    
    double get_trait_probability(unsigned person_id, 
                                 enum phased_trait pt, 
                                 unsigned locus);
    double get_transmission_probability(enum phased_trait parent);

 public :
    TraitRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2);
    virtual ~TraitRfunction() {}
    TraitRfunction(const TraitRfunction& rhs);
    TraitRfunction& operator=(const TraitRfunction& rhs);
};

#endif

