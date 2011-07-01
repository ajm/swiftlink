#ifndef LKG_LOCUSSAMPLER_H_
#define LKG_LOCUSSAMPLER_H_

using namespace std;

#include <vector>

#include "descent_graph_diff.h"
#include "descent_graph.h"
#include "peeler.h"


class Pedigree;
class Rfunction;
class SamplerRfunction;
class GeneticMap;


const unsigned int _SAMPLE_PERIOD = 100;

class LocusSampler {

    Pedigree* ped;
    GeneticMap* map;
    DescentGraph dg;
    vector<SamplerRfunction*> rfunctions;
    Peeler peel;
    
    unsigned sample_mi(unsigned allele, enum phased_trait trait, 
                       unsigned personid, unsigned locus, enum parentage parent);
    unsigned sample_homo_mi(unsigned personid, unsigned locus, enum parentage parent);
    unsigned sample_hetero_mi(unsigned allele, enum phased_trait trait);
    unsigned update_temperature(unsigned temps, unsigned current_temp);
    unsigned update_temperature_hastings(unsigned temps, unsigned current_temp);
    unsigned get_random(unsigned i);
    unsigned get_random_locus();
    double get_random();
    void init_rfunctions();
    void copy_rfunctions(const LocusSampler& rhs);
    void kill_rfunctions();
    void step(double temperature=0.0);

    
 public :
    LocusSampler(Pedigree* ped, GeneticMap* map);        
    ~LocusSampler();
    LocusSampler(const LocusSampler& rhs);
    LocusSampler& operator=(const LocusSampler& rhs);
    
    Peeler* run(unsigned iterations);
    Peeler* temper(unsigned iterations, unsigned num_temperatures);
};

#endif

