#ifndef LKG_ELOD_H_
#define LKG_ELOD_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "pedigree_parser.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "types.h"
#include "trait.h"
#include "random.h"
#include "disease_model.h"


class Elod {

    DiseaseModel dm;
    vector<Pedigree> pedigrees;
    GeneticMap map1;
    GeneticMap map2;
    struct mcmc_options options;

  public :
    Elod(const char* pedfile, struct mcmc_options opt) : 
        dm(opt.elod_frequency, opt.elod_penetrance, false), 
        pedigrees(), 
        map1(1),
        map2(1),
        options(opt) {

        
        // setup a fake map
        Snp s0("marker", 0.0);
        s0.set_minor_freq(0.5);

        Snp s1("trait", options.elod_marker_separation / 2);
        s1.set_minor_freq(options.elod_frequency);

        Snp s2("marker", options.elod_marker_separation);
        s2.set_minor_freq(0.5);

        map1.add(s0);
        map1.add(s1);
        map1.add(s2);
        map1.add_theta(map1.haldane(options.elod_marker_separation / 2));
        map1.add_theta(map1.haldane(options.elod_marker_separation / 2));

        map2.add(s0);
        map2.add(s2);
        map2.add_theta(map2.haldane(options.elod_marker_separation));

        if(! map1.sanity_check() or ! map2.sanity_check()) {
            exit(1);
        }


        // read in pedigrees
        PedigreeParser pp(pedfile, pedigrees, dm, map1);
        pp.set_ignore_genotypes();
        if(! pp.parse()) {
            exit(1);
        }

        // kill any genotype information read in
        // add UNTYPED for both loci in the fake map
        // then copy across the probabilities from the disease trait
        for(unsigned int i = 0; i < pedigrees.size(); ++i) {
            for(unsigned int j = 0; j < pedigrees[i].num_members(); ++j) {
                Person* p = pedigrees[i].get_by_index(j);

                p->clear_genotypes();

                p->add_genotype(UNTYPED);
                p->add_genotype(UNTYPED);
                p->add_genotype(UNTYPED);

                p->populate_trait_prob_cache(map1);

                p->copy_disease_probs(1);
            }
        }


        // set random seed
        init_random();
        if(options.random_filename == "") {
            seed_random_implicit();
        }
        else { 
            seed_random_explicit(options.random_filename);
        }
    }

    ~Elod() {}

    Elod(const Elod& rhs) :
        dm(rhs.dm),
        pedigrees(rhs.pedigrees),
        map1(rhs.map1),
        map2(rhs.map2),
        options(rhs.options) {}

    Elod& operator=(const Elod& rhs) {
        if(this != &rhs) {
            dm = rhs.dm;
            pedigrees = rhs.pedigrees;
            map1 = rhs.map1;
            map2 = rhs.map2;
            options = rhs.options;
        }
        return *this;
    }

    double run();
};

#endif

