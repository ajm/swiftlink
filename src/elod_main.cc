using namespace std;

#include <cstdio>
#include <cstdlib>

#include "omp.h"
#include "types.h"
#include "elod.h"


int main() {
    double frequency = 0.0001;
    double penetrance[3];

    penetrance[0] = penetrance[1] = 0.0;
    penetrance[2] = 1.0;

    struct mcmc_options opt;

    opt.peelopt_iterations = 10000; // XXX
    opt.iterations = 100000;
    opt.verbose = true;
    opt.elod_frequency = frequency;
    opt.elod_marker_separation = 0.08;

    //opt.random_filename = "random.txt";

    for(int i = 0; i < 3; ++i)
        opt.elod_penetrance[i] = penetrance[i];

    omp_set_num_threads(1);

    Elod e("ped", opt);

    printf("ELOD = %.2f\n", e.run());

    return EXIT_SUCCESS;
}

