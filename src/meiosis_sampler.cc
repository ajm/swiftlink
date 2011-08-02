using namespace std;

#include <cmath>

#include "pedigree.h"
#include "logarithms.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "descent_graph_types.h"
#include "meiosis_sampler.h"
#include "founder_allele_graph.h"


void MeiosisSampler::init_matrices() {
    matrix = new MeiosisMatrix[map->num_markers()];
}

void MeiosisSampler::copy_matrices(const MeiosisSampler& rhs) {
    copy(rhs.matrix, rhs.matrix + map->num_markers(), matrix);
}

void MeiosisSampler::kill_matrices() {
    delete[] matrix;
}

double MeiosisSampler::graph_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, unsigned value) {
    double tmp_likelihood;
    unsigned tmp = dg.get(person_id, locus, parent);
    
    dg.set(person_id, locus, parent, value);
    
    f.reset();
    
    if(f.populate(dg, locus)) {
        if(not f.likelihood(&tmp_likelihood, locus)) {        
            tmp_likelihood = LOG_ILLEGAL;
        }
    }
    else {
        tmp_likelihood = LOG_ILLEGAL;
    }
    
    dg.set(person_id, locus, parent, tmp);
    
    return tmp_likelihood;
}

void MeiosisSampler::step(DescentGraph& dg, unsigned parameter) {
    enum parentage p = get_random_meiosis();
    
    // forwards
    matrix[0][0] = graph_likelihood(dg, parameter, 0, p, 0);
    matrix[0][1] = graph_likelihood(dg, parameter, 0, p, 1);
    
    for(unsigned i = 1; i < map->num_markers(); ++i) {
        for(unsigned j = 0; j < 2; ++j) {
            // these are log likelihood, so need to be handled carefully
            matrix[i][j] = log_product( \
                                       graph_likelihood(dg, parameter, i, p, j), \
                                       log_sum( \
                                               log_product(matrix[i-1][j], map->get_theta(i-1)), \
                                               log_product(matrix[i-1][1-j], map->get_inverse_theta(i-1)) \
                                               ) \
                                       );
        }
    }
    
    // sample backwards
    // change descent graph in place
    int i = map->num_markers() - 1;
    matrix[i].normalise();
    //matrix[i].print();
    dg.set(parameter, i, p, matrix[i].sample());
    
    while(--i >= 0) {
        for(int j = 0; j < 2; ++j) {
            matrix[i][j] = log_product(matrix[i][j], ((dg.get(parameter, i+1, p) != j) ? map->get_theta(i) : map->get_inverse_theta(i)));
        }
        
        matrix[i].normalise();
        //matrix[i].print();
        dg.set(parameter, i, p, matrix[i].sample());
    }
    
    /*
    // XXX comment out when I know everything is cool
    if(not dg.likelihood()) {
        fprintf(stderr, "Error: descent graph produced by M-sampler is illegal!\n");
        abort();
    }
    */
}
