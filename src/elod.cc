using namespace std;

#include <cstdio>

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

double Elod::run() {
    double total_elod = 0.0;

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        PeelSequenceGenerator psg(&pedigrees[i], &map1, options.verbose);
        psg.build_peel_sequence(options.peelopt_iterations);

        DescentGraph dg1(&pedigrees[i], &map1);
        DescentGraph dg2(&pedigrees[i], &map2);

        LocusSampler lsampler(&pedigrees[i], &map1, &psg, 0);
        LODscores lod(&map2);
        
        Peeler peel(&pedigrees[i], &map2, &psg, &lod);
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

            //printf("%s\n", lod.debug_string().c_str());
        }

        p.finish();

        total_elod += lod.get(0,0);

        if(options.verbose)
            fprintf(stderr, "Pedigree %d ELOD = %.2f\n", i, lod.get(0,0));
    
        //for(int i = 0; i < 10; ++i) {
        //    fprintf(stderr, " %d %.2f\n", i, lod.get(0,i));
        //}
    }

    if(options.verbose)
        fprintf(stderr, "Final ELOD = %.2f\n", total_elod);

    return total_elod;
}

