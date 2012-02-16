using namespace std;

#include <cmath>

#include "types.h"
#include "pedigree.h"
#include "logarithms.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "meiosis_sampler.h"
#include "founder_allele_graph4.h"
#include "random.h"

#include "founder_allele_graph3.h"


void MeiosisSampler::reset(DescentGraph& dg) {
    #pragma omp parallel for
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        f4[i].reset(dg);
    }
}

void MeiosisSampler::find_founderallelegraph_ordering() {
    vector<bool> visited(ped->num_members(), false);
    int total = ped->num_members();
    
    // we can start by putting in all founders as there are clearly
    // no dependencies
    for(unsigned i = 0; i < ped->num_founders(); ++i) {
        seq.push_back(i);
        visited[i] = true;
        total--;
    }
    
    while(total > 0) {
        for(unsigned i = ped->num_founders(); i < ped->num_members(); ++i) {
            if(visited[i])
                continue;
        
            Person* p = ped->get_by_index(i);
            
            if(visited[p->get_maternalid()] and visited[p->get_paternalid()]) {
                seq.push_back(i);
                visited[i] = true;
                total--;
            }
        }
    }
    
    if(seq.size() != ped->num_members()) {
        fprintf(stderr, "error: founder allele sequence generation failed\n");
        abort();
    }
}

double MeiosisSampler::graph_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, unsigned value) {
    double lik4;//, lik3;
    unsigned tmp = dg.get(person_id, locus, parent);
    //bool flip = tmp != value;
    
    /*
    if(flip) {
        f4[locus].flip(person_id, parent);
    }
    */
    
    dg.set(person_id, locus, parent, value);
    
    f4[locus].reset(dg);
    
    lik4 = f4[locus].likelihood();
    
    /*
    lik3 = f3.evaluate(dg, locus);
    
    if(fabs(lik3 - lik4) > 0.000001) {
        fprintf(stderr, "%d,%d F3 %f F4 %f\n", locus, value, lik3, lik4);
        fprintf(stderr, "%s\n\n", f3.debug_string().c_str());
        abort();
    }
    */
    
    /*
    if(flip) {
        //fprintf(stderr, "2 flip\n");
        f4[locus].flip(person_id, parent);
    }
    */
    
    dg.set(person_id, locus, parent, tmp);
    
    return lik4 == 0.0 ? LOG_ZERO : log(lik4);
}

void MeiosisSampler::step(DescentGraph& dg, unsigned parameter) {
    // parameter is the founder allele
    unsigned person_id = ped->num_founders() + (parameter / 2);
    enum parentage p = static_cast<enum parentage>(parameter % 2);
    int num_markers = static_cast<int>(map->num_markers());
    int i, j;
    //int tmp, tmp2;
    
    
    #pragma omp parallel for
    for(i = 0; i < num_markers; ++i) {
        int index = i * 2;
        
        //if(parameter == 0) {
            matrix[index + 0] = graph_likelihood(dg, person_id, i, p, 0);
            matrix[index + 1] = graph_likelihood(dg, person_id, i, p, 1);
        /*
        }
        else {
            tmp = dg.get(person_id, i, p);
            tmp2 = dg.get(ped->num_founders() + ((parameter - 1) / 2), i, static_cast<enum parentage>((parameter - 1) % 2));
            
            matrix[index + tmp] = matrix[index + tmp2];
            matrix[index + (1-tmp)] = graph_likelihood(dg, person_id, i, p, 1-tmp);
        }
        */
        
        if((matrix[index] == LOG_ZERO) and (matrix[index + 1] == LOG_ZERO)) {
            fprintf(stderr, "bad: (locus = %d)\n", i);
            abort();
        }
    }
    
    
    // forwards
    for(i = 1; i < num_markers; ++i) {
        int index = i * 2;
        
        for(j = 0; j < 2; ++j) {
            /*
            matrix[index + j] = log_product( \
                                       matrix[index + j], \
                                       log_sum( \
                                               log_product(matrix[((i-1) * 2) + (1-j)], map->get_theta_log(i-1)), \
                                               log_product(matrix[((i-1) * 2) +    j ], map->get_inversetheta_log(i-1)) \
                                               ) \
                                       );
            */
            
            double prev1 = log_product(matrix[((i-1) * 2) + (1-j)], map->get_theta_log(i-1));
            double prev2 = log_product(matrix[((i-1) * 2) +    j ], map->get_inversetheta_log(i-1));
            double total = log_sum(prev1, prev2);
            prev1 -= total;
            prev2 -= total;
            
            matrix[index + j] = log_sum(
                                    log_product(matrix[index + j], prev1), 
                                    log_product(matrix[index + j], prev2)
                                );
        }
    }
    
    // sample backwards
    // change descent graph in place
    i = num_markers - 1;
    dg.set(person_id, i, p, sample(i));
    /*
    tmp = sample(i);
    tmp2 = dg.get(person_id, i, p);
    
    if(tmp != tmp2) {
        dg.set(person_id, i, p, tmp);
        //fprintf(stderr, "3 flip (%d)\n", i);
        f4[i].flip(person_id, p);
    }
    */
    
    while(--i >= 0) {
        int index = i * 2;
        
        for(j = 0; j < 2; ++j) {
            double next = (dg.get(person_id, i+1, p) != j) ? map->get_theta_log(i) : map->get_inversetheta_log(i);
        
            matrix[index + j] = log_product(matrix[index + j], next);
        }
        
        /*
        tmp = sample(i);
        tmp2 = dg.get(person_id, i, p);
        
        if(tmp2 != tmp) {
            dg.set(person_id, i, p, tmp);
            //fprintf(stderr, "4 flip (%d)\n", i);
            //f4[i].flip(person_id, p);
        }
        */
        
        dg.set(person_id, i, p, sample(i));
    }
}

unsigned MeiosisSampler::sample(int locus) {
    int index = locus * 2;
    
    if((matrix[index] == LOG_ZERO) and (matrix[index + 1] == LOG_ZERO)) {
        fprintf(stderr, "error: m-sampler probability distribution sums to zero (locus = %d)\n", locus);
        abort();
    }
    
    if(matrix[index] == LOG_ZERO) {
        return 1;
    }
    
    if(matrix[index+1] == LOG_ZERO) {
        return 0;
    }
    
    double r = get_random();
    
    return (log(r) < (matrix[index] - log_sum(matrix[index], matrix[index + 1]))) ? 0 : 1;
}

