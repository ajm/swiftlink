using namespace std;

#include <cmath>

#include "types.h"
#include "pedigree.h"
#include "logarithms.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "meiosis_sampler.h"
#include "founder_allele_graph2.h"
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
    /*
    f.set_locus(locus);
    f.reset();
    
    tmp_likelihood = f.populate(dg) ? f.likelihood() : LOG_ILLEGAL;
    */
    
    fag.reset();
    if(fag.populate(dg, locus)) {
        //fag.print();
        if(not fag.likelihood(&tmp_likelihood, locus)) {        
            tmp_likelihood = LOG_ILLEGAL;
        }
    }
    else {
        tmp_likelihood = LOG_ILLEGAL;
    }
        
    dg.set(person_id, locus, parent, tmp);
    
    /*
    if(tmp_likelihood != tmp_likelihood2) {
        fprintf(stderr, "old = %e, new = %e, illegal = %e\n", tmp_likelihood2, tmp_likelihood, LOG_ILLEGAL);
        fag.print();
        fprintf(stderr, "%s", f.debug_string().c_str());
        abort();
    }
    */
    
    return tmp_likelihood;
}

/*
double MeiosisSampler::initial_likelihood(DescentGraph& dg, unsigned locus) {

    f.set_locus(locus);
    f.reset();
    
    return f.populate(dg) ? f.likelihood() : LOG_ILLEGAL;
}

void MeiosisSampler::incremental_likelihood(DescentGraph& dg, 
                                            unsigned person_id, unsigned locus, enum parentage parent, 
                                            double* meiosis0, double* meiosis1) {
    
    int current = dg.get(person_id, locus, parent);
    double tmp1;
    double tmp2;
    double tmp3 = 0.0;
    
    fprintf(stderr, "incremental_likelihood start\n");
    
    fprintf(stderr, "tmp1\n");
    tmp1 = initial_likelihood(dg, locus);
    
    if(tmp1 == LOG_ILLEGAL) {
        dg.flip_bit(person_id, locus, parent);
        fprintf(stderr, "tmp2 (1)\n");
        tmp2 = initial_likelihood(dg, locus);
    }
    else {
        fprintf(stderr, "tmp2 (2)\n");
        tmp2 = f.reevaluate(dg, person_id, locus, parent);
        
        fprintf(stderr, "tmp3\n");
        tmp3 = initial_likelihood(dg, locus);
        
        fprintf(stderr, "NO INCREMENTAL:\n%s\n", f.debug_string().c_str());
        
        if(tmp2 != tmp3) {
            fprintf(stderr, "reevaluate() = %e, actual = %e, prev = %e (%d --> %d)\n", tmp2, tmp3, tmp1, current, dg.get(person_id, locus, parent));
            fprintf(stderr, "%s\n", f.debug_string().c_str());
            abort();
        }
    }
    
    *meiosis0 = (current == 0) ? tmp1 : tmp2;
    *meiosis1 = (current == 0) ? tmp2 : tmp1;
    
    fprintf(stderr, "incremental_likelihood end (tmp1 = %e, tmp2 = %e tmp3 = %e)\n\n", tmp1, tmp2, tmp3);
}
*/

void MeiosisSampler::step(DescentGraph& dg, unsigned parameter) {
    //enum parentage p = get_random_meiosis();
    //double meiosis0, meiosis1;
    
    // parameter is the founder allele
    unsigned person_id = ped->num_founders() + (parameter / 2);
    enum parentage p = static_cast<enum parentage>(parameter % 2);
    
    //fprintf(stderr, "%d %s\n", int(person_id), p == MATERNAL ? "M" : "P");
    
    // forwards
    matrix[0][0] = graph_likelihood(dg, person_id, 0, p, 0);
    matrix[0][1] = graph_likelihood(dg, person_id, 0, p, 1);
    //incremental_likelihood(dg, parameter, 0, p, &meiosis0, &meiosis1);
    //matrix[0][0] = meiosis0;
    //matrix[0][1] = meiosis1;
    
    for(unsigned i = 1; i < map->num_markers(); ++i) {
        for(unsigned j = 0; j < 2; ++j) {
            // these are log likelihood, so need to be handled carefully
            matrix[i][j] = log_product( \
                                       graph_likelihood(dg, person_id, i, p, j), \
                                       log_sum( \
                                               log_product(matrix[i-1][j],   map->get_theta_log(i-1)), \
                                               log_product(matrix[i-1][1-j], map->get_inversetheta_log(i-1)) \
                                               ) \
                                       );
        }
        
        //double theta = map->get_theta(i-1);
        //double antitheta = map->get_inverse_theta(i-1);
        
        //incremental_likelihood(dg, parameter, i, p, &meiosis0, &meiosis1);
        
        //matrix[i][0] = log_product(meiosis0, log_sum(log_product(matrix[i-1][0], theta), log_product(matrix[i-1][1], antitheta)));
        //matrix[i][1] = log_product(meiosis1, log_sum(log_product(matrix[i-1][1], theta), log_product(matrix[i-1][0], antitheta)));
    }
    
    // sample backwards
    // change descent graph in place
    int i = map->num_markers() - 1;
    matrix[i].normalise();
    //matrix[i].print();
    dg.set(person_id, i, p, matrix[i].sample());
    
    while(--i >= 0) {
        for(int j = 0; j < 2; ++j) {
            matrix[i][j] = log_product(matrix[i][j], ((dg.get(person_id, i+1, p) != j) ? map->get_theta_log(i) : map->get_inversetheta_log(i)));
        }
        
        matrix[i].normalise();
        //matrix[i].print();
        dg.set(person_id, i, p, matrix[i].sample());
    }
    
    /*
    // XXX comment out when I know everything is cool
    if(not dg.likelihood()) {
        fprintf(stderr, "Error: descent graph produced by M-sampler is illegal!\n");
        abort();
    }
    */
}

