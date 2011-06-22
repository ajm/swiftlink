#ifndef LKG_SAMPLERRFUNCTION_H_
#define LKG_SAMPLERRFUNCTION_H_

#include <vector>

#include "trait.h"
#include "rfunction.h"

class DescentGraph;
class Pedigree;
class GeneticMap;

class SamplerRfunction : public Rfunction {
    
    double get_recombination_probability(DescentGraph* dg, 
                                         unsigned locus,
                                         unsigned person_id, 
                                         int maternal_allele, 
                                         int paternal_allele);
    
    double get_trait_probability(unsigned person_id, 
                                 enum phased_trait pt, 
                                 unsigned locus);

 public :
    SamplerRfunction(PeelOperation po, 
                     Pedigree* p, 
                     GeneticMap* m, 
                     vector<Rfunction*>& previous_functions, 
                     unsigned index);
    virtual ~SamplerRfunction() {}
    SamplerRfunction(const Rfunction& rhs);
    SamplerRfunction& operator=(const Rfunction& rhs);
};

#endif
