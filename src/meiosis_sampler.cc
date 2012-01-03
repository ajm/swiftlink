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
#include "founder_allele_graph4.h"


void MeiosisSampler::init_matrices() {
    matrix = new MeiosisMatrix[map->num_markers()];
}

void MeiosisSampler::copy_matrices(const MeiosisSampler& rhs) {
    copy(rhs.matrix, rhs.matrix + map->num_markers(), matrix);
}

void MeiosisSampler::kill_matrices() {
    delete[] matrix;
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
    //double tmp_likelihood;
    unsigned tmp = dg.get(person_id, locus, parent);
    
    dg.set(person_id, locus, parent, value);
        
    // FounderAlleleGraph3
    //tmp_likelihood = f3.evaluate(dg, locus); // XXX
    //printf("%s\n", f3.debug_string().c_str());
    
    double tmp_fag4 = f4.init_likelihood(dg, locus);
    tmp_fag4 = tmp_fag4 == 0.0 ? LOG_ILLEGAL : log(tmp_fag4);
    
    dg.set(person_id, locus, parent, tmp);
    /*
    if(tmp_likelihood != tmp_fag4) {
        fprintf(stderr, "old = %e, new = %e, illegal = %e\n", tmp_likelihood, tmp_fag4, LOG_ILLEGAL);
        fprintf(stderr, "%s\n", f3.debug_string().c_str());
        abort();
    }
    */
    return tmp_fag4; //tmp_likelihood;
}

void MeiosisSampler::step(DescentGraph& dg, unsigned parameter) {
    // parameter is the founder allele
    unsigned person_id = ped->num_founders() + (parameter / 2);
    enum parentage p = static_cast<enum parentage>(parameter % 2);
    
    //fprintf(stderr, "%d %s\n", int(person_id), p == MATERNAL ? "M" : "P");
    
    // forwards
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        if(parameter == 0) {
            matrix[i][0] = graph_likelihood(dg, person_id, i, p, 0);
            matrix[i][1] = graph_likelihood(dg, person_id, i, p, 1);
        }
        else {
            int tmp = dg.get(person_id, i, p);
            int tmp2 = dg.get(ped->num_founders() + ((parameter - 1) / 2), i, static_cast<enum parentage>((parameter - 1) % 2));
            matrix[i][tmp] = matrix[i][tmp2];
            matrix[i][1-tmp] = graph_likelihood(dg, person_id, i, p, 1-tmp);
        }
    }
    
    for(unsigned i = 1; i < map->num_markers(); ++i) {
        for(unsigned j = 0; j < 2; ++j) {
            // these are log likelihood, so need to be handled carefully
            matrix[i][j] = log_product( \
                                       matrix[i][j], /*graph_likelihood(dg, person_id, i, p, j),*/ \
                                       log_sum( \
                                               log_product(matrix[i-1][j],   map->get_theta_log(i-1)), \
                                               log_product(matrix[i-1][1-j], map->get_inversetheta_log(i-1)) \
                                               ) \
                                       );
        }
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

