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


void MeiosisSampler::reset(DescentGraph& dg, unsigned int parameter) {
    unsigned person_id = ped->num_founders() + (parameter / 2);
    enum parentage p = static_cast<enum parentage>(parameter % 2);
    
    #pragma omp parallel for
    for(unsigned int i = 0; i < map->num_markers(); ++i) {
        int index = i * 2;
        int meiosis = dg.get(person_id, i, p);
        
        matrix[index + meiosis] = graph_likelihood(dg, person_id, i, p, meiosis);
        
        if(matrix[index + meiosis] == 0.0) {
            fprintf(stderr, "error: illegal descent graph given to m-sampler (%s:%d)\n", __FILE__, __LINE__);
            abort();
        }
    }
    
    last_parameter = parameter;
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
        
    dg.set(person_id, locus, parent, value);
    
    f4[locus].reset(dg);
    lik4 = f4[locus].likelihood();
    
    dg.set(person_id, locus, parent, tmp);
    
    return lik4;
}

void MeiosisSampler::step(DescentGraph& dg, unsigned int parameter) {
    // parameter is the founder allele
    unsigned person_id = ped->num_founders() + (parameter / 2);
    enum parentage p = static_cast<enum parentage>(parameter % 2);
    
    unsigned last_id = ped->num_founders() + (last_parameter / 2);
    enum parentage last_p = static_cast<enum parentage>(last_parameter % 2);
    
    int num_markers = static_cast<int>(map->num_markers());
    
    
    #pragma omp parallel for
    for(int i = 0; i < num_markers; ++i) {
        int index = i * 2;
        
        //matrix[index + 0] = graph_likelihood(dg, person_id, i, p, 0);
        //matrix[index + 1] = graph_likelihood(dg, person_id, i, p, 1);

        
        int tmp = dg.get(person_id, i, p);
        int tmp2 = dg.get(last_id, i, last_p);
        
        matrix[index + tmp] = matrix[index + tmp2];
        matrix[index + (1-tmp)] = graph_likelihood(dg, person_id, i, p, 1-tmp);
        
        
        if((matrix[index] == 0.0) and (matrix[index+1] == 0.0)) {
            fprintf(stderr, "error: illegal descent graph given to m-sampler (%s:%d)\n", __FILE__, __LINE__);
            abort();
        }
    }
    
    
    double total = matrix[0] + matrix[1];
    
    matrix[0] /= total;
    matrix[1] /= total;
    
    // forwards
    for(int i = 1; i < num_markers; ++i) {
        int index = i * 2;
        
        for(int j = 0; j < 2; ++j) {
            matrix[index + j] *= \
                    ((matrix[((i-1) * 2) + (1-j)] * map->get_theta(i-1)) + \
                     (matrix[((i-1) * 2) +    j ] * map->get_inversetheta(i-1)));
        }
        
        double total = matrix[index] + matrix[index + 1];
        
        matrix[index] /= total;
        matrix[index + 1] /= total;
    }
    
    // sample backwards
    // change descent graph in place
    int i = num_markers - 1;
    dg.set(person_id, i, p, sample(i));
    
    while(--i >= 0) {
        int index = i * 2;
        
        for(int j = 0; j < 2; ++j) {
            double next = (dg.get(person_id, i+1, p) != j) ? map->get_theta(i) : map->get_inversetheta(i);
            matrix[index + j] *= next;
        }
        
        dg.set(person_id, i, p, sample(i));
    }
    
    last_parameter = parameter;
}

unsigned MeiosisSampler::sample(int locus) {
    int index = locus * 2;
    
    if(matrix[index] == 0.0)
        return 1;
        
    if(matrix[index + 1] == 0.0)
        return 0;
        
    return (get_random() < (matrix[index] / (matrix[index] + matrix[index + 1]))) ? 0 : 1;
}

