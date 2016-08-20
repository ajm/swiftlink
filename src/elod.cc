#include <cstdio>
#include <vector>
#include <numeric>

#include "elod.h"
#include "lod_score.h"
#include "peel_sequence_generator.h"
#include "peeler.h"
#include "locus_sampler2.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "progress.h"
#include "types.h"

using namespace std;


double Elod::run() {
    vector<double> elods;

    fprintf(stderr, "\nELOD parameters:\n"
                    "\tpenetrance = %.2f:%.2f:%.2f\n"
                    "\tseparation = %.2f\n"
                    "\ttrait freq = %.2e\n"
                    "\treplicates = %d\n"
                    "\tsex-linked = %s\n\n", 
                    options.elod_penetrance[0],
                    options.elod_penetrance[1],
                    options.elod_penetrance[2],
                    options.elod_marker_separation,
                    options.elod_frequency,
                    options.elod_replicates,
                    options.sex_linked ? "true" : "false");

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        PeelSequenceGenerator psg(&pedigrees[i], &map1, options.sex_linked, options.verbose);
        psg.build_peel_sequence(options.peelopt_iterations);

        DescentGraph dg1(&pedigrees[i], &map1, options.sex_linked);
        DescentGraph dg2(&pedigrees[i], &map2, options.sex_linked);

        LocusSampler lsampler(&pedigrees[i], &map1, &psg, 0, options.sex_linked);
        LODscores lod(&map2);
        
        Peeler peel(&pedigrees[i], &map2, &psg, &lod, options.sex_linked);
        peel.calc_trait_prob();
        peel.set_locus(0);

        lod.set_trait_prob(peel.calc_trait_prob());


        Progress p("ELOD: ", options.elod_replicates);

        for(int j = 0; j < options.elod_replicates; ++j) {
            lsampler.start_from(dg1, 1);

            dg2.copy_locus(dg1, 0, 0);
            dg2.copy_locus(dg1, 2, 1);

            peel.process(&dg2);
        
            p.increment();
        }

        p.finish();

        elods.push_back(lod.get(0,0));
    }


    double total_elod = 0.0;
    fprintf(stderr, "\n%10s |%12s\n", "Pedigree", "ELOD");
    fprintf(stderr, "-----------|------------\n");

    for(unsigned int i = 0; i < elods.size(); ++i) {
        fprintf(stderr, "%10d |%12f\n", i, elods[i]);
        total_elod += elods[i];
    }

    fprintf(stderr, "%10s |%12f\n", "Total", total_elod);

    return total_elod;
}

