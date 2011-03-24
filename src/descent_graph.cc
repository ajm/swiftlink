using namespace std;

#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "descent_graph.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "genotype.h"
#include "founder_allele_graph.h"
#include "elimination.h"


DescentGraph::DescentGraph(Pedigree* ped, GeneticMap* map) 
    : ped(ped), map(map), prob(0.0) {
	graph_size = 2 * ped->num_members();

    marker_transmission = log(0.5) * 
		(ped->num_markers() * 2 * (ped->num_members() - ped->num_founders()));
	
    data = new char[graph_size * ped->num_markers()];
}

DescentGraph::DescentGraph(const DescentGraph& d) 
	: ped(d.ped), map(d.map), prob(d.prob), 
	  marker_transmission(d.marker_transmission),
      graph_size(d.graph_size) {

    unsigned int data_length = graph_size * ped->num_markers();
	data = new char[data_length];
    copy(d.data, d.data + data_length, data);
}

DescentGraph::~DescentGraph() {
	delete[] data;
}

DescentGraph& DescentGraph::operator=(const DescentGraph& d) {
    
	if(&d != this) {
		ped = d.ped;
		map = d.map;
		prob = d.prob;
		graph_size = d.graph_size;
		marker_transmission = d.marker_transmission;
        
        delete[] data;

        unsigned int data_length = graph_size * ped->num_markers();
    
    	data = new char[data_length];
        copy(d.data, d.data + data_length, data);
	}

	return *this;
}

void DescentGraph::copy_from(DescentGraph& d, unsigned start, unsigned end) {
    //memcpy(data + start, d.data + start, end - start);
    copy(d.data + start, d.data + end, data);
}

unsigned DescentGraph::data_length() {
    return graph_size * ped->num_markers();
}

void DescentGraph::flip_bit(unsigned i) {
    data[i] = (data[i] == 0) ? 1 : 0; 
}

void DescentGraph::flip_bit(unsigned person_id, unsigned locus, enum parentage p) {
	flip_bit(_offset(person_id, locus, p));
}

char DescentGraph::get_bit(unsigned i) {
    return data[i];
}

void DescentGraph::set_bit(unsigned i, char b) {
    data[i] = b;
}

bool DescentGraph::random_descentgraph() {
    GenotypeElimination ge(ped);
    return ge.random_descentgraph(*this);
}

int DescentGraph::_offset(unsigned person_id, unsigned locus, enum parentage p) {
	return (graph_size * locus) + (person_id * 2) + p;
}

char DescentGraph::get(unsigned person_id, unsigned locus, enum parentage p) {
    return data[_offset(person_id, locus, p)];
}

void DescentGraph::set(unsigned person_id, unsigned locus, enum parentage p, char value) {
    data[_offset(person_id, locus, p)] = value;
}

// this only works because the founders are guarenteed to be at the start of
// the list of family members
int DescentGraph::_founder_allele(unsigned person_id, enum parentage p) {
    return (person_id * 2) + p;
}

int DescentGraph::get_founderallele(unsigned person_id, unsigned locus, enum parentage p) {
    unsigned current = person_id;
    enum parentage parent_allele = p;
	Person* per;
    
    while(1) {
		per = ped->get_by_index(current);

        if(per->isfounder()) {
            return _founder_allele(current, parent_allele);
        }
        
        if( parent_allele == PATERNAL ) {
            parent_allele = static_cast<enum parentage>(get(current, locus, parent_allele));
            current = per->get_paternalid();
        }
        else {
            parent_allele = static_cast<enum parentage>(get(current, locus, parent_allele));
            current = per->get_maternalid();
		}
    }
}

// transmission prob can be low, but never zero
// sum of prior prob can be zero
bool DescentGraph::likelihood(double *p) {
    double trans, prior;
	
    trans = _transmission_prob();
	
    if(not _sum_prior_prob(&prior)) {
        // good idea?
        // this is technically the same as a -INFINITY constant
        // but more portable, C requires some C99 #define, but I 
        // guess people just use DBL_MIN
        *p = prob = LOG_ILLEGAL;
        return false;
	}
    
//    printf("trans = %f\nprior = %f\n", trans / log(10), prior / log(10));
	
    *p = prob = trans + prior;
	
    return true;
}

bool DescentGraph::likelihood() {
	double p;

	return likelihood(&p);
}

bool DescentGraph::haplotype_likelihood(double *p) {
    double trans, prior;
	
    trans = _transmission_prob();
	
    if(not _best_prior_prob(&prior)) {
        *p = prob = LOG_ILLEGAL;
        return false;
	}
	
    *p = prob = trans + prior;

	//printf("trans = %.4f, prior = %.4f\n", trans / log(10), prior / log(10));
	
    return true;
}

bool DescentGraph::haplotype_likelihood() {
	double p;
    
	return haplotype_likelihood(&p);
}

double DescentGraph::_transmission_prob() {
    return marker_transmission + _recombination_prob();
}

// NOTE: descent graph storage not very cache friendly?
double DescentGraph::_recombination_prob() {
    double tmp = 0.0;
    int last, current;
	Person* p;
	
    // i = member of family
    // j = mother and father of i
    // k = current loci being considered
	
    for(unsigned i = 0; i < ped->num_members(); ++i) { // every person
		p = ped->get_by_index(i);

        if( p->isfounder() )
            continue;
        
        for(unsigned j = 0; j < 2; ++j ) { // mother and father
            last = get(i, 0, static_cast<enum parentage>(j));
            
            for(unsigned k = 1; k < ped->num_markers(); ++k ) { // every loci
                current = get(i, k, static_cast<enum parentage>(j));
				
                tmp += (last != current) ? \
						map->get_theta(k-1) : 
						map->get_inverse_theta(k-1) ;
                
                last = current;
            }
        }
    }
    
    return tmp;
}

// TODO: refactor linkage/haplotype differences
bool DescentGraph::_sum_prior_prob(double *prob) {
    double tmp_prob;
	double return_prob = 0.0;
	FounderAlleleGraph fag(map, ped);
    
	for(unsigned i = 0; i < ped->num_markers(); ++i) {
		fag.reset(); // strictly necessary?
        
        if(not fag.populate(*this, i)) {
			//printf("bad fag populate\n");
            return false;
        }
		
		if(not fag.likelihood(&tmp_prob, i)) {
            //printf("bad fag likelihood\n");
			return false;
        }
		
		return_prob += tmp_prob;
    }
	
	*prob = return_prob;    
	return true;
}

bool DescentGraph::_best_prior_prob(double *prob) {
    double tmp_prob;
	double return_prob = 0.0;
	FounderAlleleGraph fag(map, ped);
    
	for(unsigned i = 0; i < ped->num_markers(); ++i) {
		fag.reset(); // strictly necessary?
        
        if(not fag.populate(*this, i)) {
			//printf("bad fag populate\n");
            return false;
        }
		
		if(not fag.likelihood(&tmp_prob, i)) {
            //printf("bad fag likelihood\n");
			return false;
        }
		
		return_prob += fag.descentstate_likelihood(i);
    }
	
	*prob = return_prob;    
	return true;
}

void DescentGraph::print() {
	int pat, mat;
	Person* p;
	
    //printf("DescentGraph: ");
    fprintf(stderr, "DescentGraph: ");
    for(int locus = 0; locus < int(ped->num_markers()); ++locus) {
        for(unsigned i = 0; i  < ped->num_members(); ++i) {
            p = ped->get_by_index(i);
                        
            if(not p->isfounder()) {
                mat = get(i, locus, MATERNAL);
                pat = get(i, locus, PATERNAL);
                
                //printf("%d%d", mat, pat);
                fprintf(stderr, "%d%d", mat, pat);
            }
        }
    }
    //printf("\n");
    fprintf(stderr, "\n");
}

double DescentGraph::trans_prob() {
    return -4.1600831433636838; // _transmission_prob();
}

