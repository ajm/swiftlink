using namespace std;

#include <cmath>

#include "types.h"
#include "pedigree.h"
#include "logarithms.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "meiosis_sampler.h"
#include "founder_allele_graph4.h"


void MeiosisSampler::reset(DescentGraph& dg) {
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
    double lik4;
    unsigned tmp = dg.get(person_id, locus, parent);
    bool flip = tmp != value;
    
    if(flip)
        f4[locus].flip(person_id, parent);
    
    lik4 = f4[locus].likelihood();
    
    if(flip)
        f4[locus].flip(person_id, parent);
    
    return lik4 == 0 ? LOG_ZERO : log(lik4);
}

void MeiosisSampler::step(DescentGraph& dg, unsigned parameter) {
    // parameter is the founder allele
    unsigned person_id = ped->num_founders() + (parameter / 2);
    enum parentage p = static_cast<enum parentage>(parameter % 2);
    int num_markers = static_cast<int>(map->num_markers());
    int i, j;
    int tmp, tmp2;
    
    // forwards
    #pragma omp parallel for
    for(i = 0; i < num_markers; ++i) {
        int index = i * 2;
        
        if(parameter == 0) {
            matrix[index + 0] = graph_likelihood(dg, person_id, i, p, 0);
            matrix[index + 1] = graph_likelihood(dg, person_id, i, p, 1);
        }
        else {
            tmp = dg.get(person_id, i, p);
            tmp2 = dg.get(ped->num_founders() + ((parameter - 1) / 2), i, static_cast<enum parentage>((parameter - 1) % 2));
            
            matrix[index + tmp] = matrix[index + tmp2];
            matrix[index + (1-tmp)] = graph_likelihood(dg, person_id, i, p, 1-tmp);
        }
    }
    
    for(i = 1; i < num_markers; ++i) {
        for(j = 0; j < 2; ++j) {
            /*
            matrix[(i * 2) + j] *= \
                                ( \
                                    (matrix[((i-1) * 2) + j]     * map->get_theta(i-1)) + \
                                    (matrix[((i-1) * 2) + (1-j)] * map->get_inversetheta(i-1)) \
                                );
            */
            matrix[(i * 2) + j] = log_product( \
                                       matrix[(i * 2) + j], \
                                       log_sum( \
                                               log_product(matrix[((i-1) * 2) + j],     map->get_theta_log(i-1)), \
                                               log_product(matrix[((i-1) * 2) + (1-j)], map->get_inversetheta_log(i-1)) \
                                               ) \
                                       );
        }
    }
    
    // sample backwards
    // change descent graph in place
    i = num_markers - 1;
    tmp = sample(i);
    tmp2 = dg.get(person_id, i, p);
    
    if(tmp != tmp2) {
        dg.set(person_id, i, p, tmp);
        f4[i].flip(person_id, p);
    }
    
    while(--i >= 0) {
        int index = i * 2;
        for(j = 0; j < 2; ++j) {
            /*
            matrix[index + j] *= ((dg.get(person_id, i+1, p) != j) ? map->get_theta(i) : map->get_inversetheta(i));
            */
            matrix[index + j] = log_product(matrix[index + j], ((dg.get(person_id, i+1, p) != j) ? map->get_theta_log(i) : map->get_inversetheta_log(i)));
        }
        
        tmp = sample(i);
        tmp2 = dg.get(person_id, i, p);
        
        if(tmp2 != tmp) {
            dg.set(person_id, i, p, tmp);
            f4[i].flip(person_id, p);
        }
    }
}

unsigned MeiosisSampler::sample(int locus) {
    int index = locus * 2;
    
    if((matrix[index] == LOG_ZERO) and (matrix[index + 1] == LOG_ZERO)) {
        fprintf(stderr, "error: m-sampler probability distribution sums to zero\n");
        abort();
    }
    
    if(matrix[index] == LOG_ZERO) {
        return 1;
    }
        
    if(matrix[index+1] == LOG_ZERO) {
        return 0;
    }
    
    return (log(get_random()) < (matrix[index] - log_sum(matrix[index], matrix[index + 1]))) ? 0 : 1;
}

