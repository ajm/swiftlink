#ifndef LKG_MARKOVCHAIN_H_
#define LKG_MARKOVCHAIN_H_

using namespace std;

#include <cstdio>

#include "types.h"
#include "genetic_map.h"
#include "meiosis_sampler.h"
#include "locus_sampler2.h"
#include "peeler.h"

class Pedigree;
class PeelSequenceGenerator;
class DescentGraph;
class LODscores;
#ifdef USE_CUDA
class GPULodscores;
#endif

class MarkovChain {
    
    Pedigree* ped;
    GeneticMap map;
    PeelSequenceGenerator* psg;
    struct mcmc_options options;
    
    LODscores* lod;
#ifdef USE_CUDA
    GPULodscores* gpulod;
#endif
    vector<Peeler*> peelers;
    vector<LocusSampler*> lsamplers;
    MeiosisSampler msampler;
    vector<int> l_ordering;
    vector<int> m_ordering;

    FILE* coda_filehandle;
    int seq_num;

    double temperature;

    void _init();
    void _kill();

 public :
    MarkovChain(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, struct mcmc_options options, int sequence_num, double temp=1.0) :
        ped(ped), 
        map(*map), 
        psg(psg),
        options(options),
        lod(0),
#ifdef USE_CUDA
        gpulod(0),
#endif
        peelers(),
        lsamplers(),
        msampler(ped, map, options.sex_linked),
        l_ordering(),
        m_ordering(),
        coda_filehandle(NULL),
        seq_num(sequence_num),
        temperature(temp) {
    
        _init();
    }
    
    ~MarkovChain() {
        _kill();
    }
    
    MarkovChain(const MarkovChain& rhs) :
        ped(rhs.ped), 
        map(rhs.map),
        psg(rhs.psg), 
        options(rhs.options),
        lod(rhs.lod),
#ifdef USE_CUDA
        gpulod(rhs.gpulod),
#endif
        peelers(rhs.peelers),
        lsamplers(rhs.lsamplers),
        msampler(rhs.msampler), 
        l_ordering(rhs.l_ordering),
        m_ordering(rhs.m_ordering),
        temperature(rhs.temperature) {}
    
    MarkovChain& operator=(const MarkovChain& rhs) {
        if(this != &rhs) {
            ped = rhs.ped;
            map = rhs.map;
            psg = rhs.psg;
            options = rhs.options;
            lod = rhs.lod;
#ifdef USE_CUDA
            gpulod = rhs.gpulod;
#endif
            peelers = rhs.peelers;
            lsamplers = rhs.lsamplers;
            msampler = rhs.msampler;
            l_ordering = rhs.l_ordering;
            m_ordering = rhs.m_ordering;
            temperature = rhs.temperature;
        }
        return *this;
    }
 
    void step(DescentGraph& dg, int start_iteration, int step_size);
    LODscores* get_result() {
        if(temperature != 1.0) {
            fprintf(stderr, "error: only the coldest chain can be used!\n");
            abort();
        }
        return lod;
    }
    double get_likelihood(DescentGraph& dg) {
        return dg.get_likelihood2(&map);
    }

    LODscores* run(DescentGraph& dg);
    bool noninterferring(vector<int>& x, int val);
};

#endif

