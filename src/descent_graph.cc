using namespace std;

#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

#include "types.h"
#include "descent_graph.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "founder_allele_graph.h"
#include "elimination.h"
#include "founder_allele_graph2.h"


DescentGraph::DescentGraph(Pedigree* ped, GeneticMap* map) : 
    data(NULL),
    ped(ped), 
    map(map), 
    prob(0.0), 
    marker_transmission(log(0.5) * (2 * (ped->num_members() - ped->num_founders()))), 
    graph_size(2 * ped->num_members()),
    recombinations(-1) {
    
    data = new int[graph_size * map->num_markers()];
    
    for(unsigned i = 0; i < (graph_size * map->num_markers()); ++i) {
        data[i] = 0;
    }    
}

DescentGraph::DescentGraph(const DescentGraph& d) : 
    data(NULL),
    ped(d.ped), 
    map(d.map), 
    prob(d.prob), 
	marker_transmission(d.marker_transmission),
    graph_size(d.graph_size),
    recombinations(d.recombinations) {

    unsigned int data_length = graph_size * map->num_markers();
    
	data = new int[data_length];
    copy(d.data, 
         d.data + data_length, 
         data);
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
                
        copy(d.data, 
             d.data + (graph_size * map->num_markers()), 
             data);
    }

	return *this;
}

bool DescentGraph::random_descentgraph() {
    GenotypeElimination ge(ped);
    return ge.random_descentgraph(*this);
}

/*
void DescentGraph::copy_from(DescentGraph& d, unsigned start, unsigned end) {
    copy(d.data + start, d.data + end, data);
}
*/

int DescentGraph::get_bit(unsigned i) const {
    return data[i];
}

void DescentGraph::set_bit(unsigned i, int b) {
    data[i] = b;
}

void DescentGraph::flip_bit(unsigned i) {
    data[i] = (data[i] == 0) ? 1 : 0; 
}

/*
int DescentGraph::_offset(unsigned person_id, unsigned locus, enum parentage p) const {
	return (graph_size * locus) + (person_id * 2) + p;
}

int DescentGraph::get(unsigned person_id, unsigned locus, enum parentage p) const {
    return data[_offset(person_id, locus, p)];
}
*/

void DescentGraph::set(unsigned person_id, unsigned locus, enum parentage p, int value) {
    data[_offset(person_id, locus, p)] = value;
}

void DescentGraph::flip(unsigned person_id, unsigned locus, enum parentage p) {
	flip_bit(_offset(person_id, locus, p));
}

// this only works because the founders are guarenteed to be at the start of
// the list of family members
int DescentGraph::_founder_allele(unsigned person_id, enum parentage p) const {
    return (person_id * 2) + p;
}

int DescentGraph::get_founderallele(unsigned person_id, unsigned locus, enum parentage p) const {
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

double DescentGraph::get_likelihood() {
    prob = log_product(_transmission_prob(), _sum_prior_prob());
    return prob;
}

double DescentGraph::get_haplotype_likelihood() {
    prob = log_product(_transmission_prob(), _best_prior_prob());
    return prob;
}

double DescentGraph::_transmission_prob() {
    return marker_transmission + _recombination_prob();
}

double DescentGraph::_recombination_prob() {
    double tmp = 0.0;
    
    recombinations = 0;
	
	for(unsigned i = 0; i < (map->num_markers() - 1); ++i) { // every loci	
        tmp += get_recombination_prob(i, true);
    }
    
    return tmp;
}

double DescentGraph::get_recombination_prob(unsigned locus, bool count_crossovers) {
    double tmp = 0.0;
	enum parentage parent;
	bool crossover;
	Person* p;
	
    double theta = map->get_theta_log(locus);
    double antitheta = map->get_inversetheta_log(locus);
	    
    for(unsigned i = 0; i < ped->num_members(); ++i) { // every person
        p = ped->get_by_index(i);
        
        if(p->isfounder())
            continue;
        
        for(unsigned j = 0; j < 2; ++j) { // mother and father
            parent = static_cast<enum parentage>(j);
            crossover = get(i, locus, parent) != get(i, locus + 1, parent);
            
            if(count_crossovers and crossover) {
                ++recombinations;
            }
				
            tmp += crossover ? theta : antitheta ;
        }
    }
    
    return tmp;
}

double DescentGraph::_sum_prior_prob() {
    /*
    double tmp_prob;
	double return_prob = 0.0;
	FounderAlleleGraph fag(ped, map);
    
	for(unsigned i = 0; i < map->num_markers(); ++i) {
		fag.reset();
        
        if(not fag.populate(*this, i)) {
            fprintf(stderr, "error: could not populate locus %d\n", int(i));
            fag.print();
            return LOG_ZERO;
        }
		
		if(not fag.likelihood(&tmp_prob, i)) {
            fprintf(stderr, "error: could not likelihood locus %d\n", int(i));
			fag.print();
            return LOG_ZERO;
        }
        		
		return_prob += tmp_prob;
    }
    
    return return_prob;
    */
    
    
    double tmp_prob;
	double return_prob = 0.0;
    FounderAlleleGraph2 f(ped, map, 0);
    
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        f.set_locus(i);
        f.reset();
        
        if(not f.populate(*this)) {
            fprintf(stderr, "error: could not populate locus %d\n", int(i));
            fprintf(stderr, "%s\n", f.debug_string().c_str());
            return LOG_ZERO;
        }
		
		if((tmp_prob = f.likelihood()) == LOG_ZERO) {
            fprintf(stderr, "error: could not likelihood locus %d\n", int(i));
            fprintf(stderr, "%s\n", f.debug_string().c_str());
			return LOG_ZERO;
        }
        
		return_prob += tmp_prob;
    }
    
    return return_prob;
    
}

double DescentGraph::_best_prior_prob() {
    double tmp_prob;
	double return_prob = 0.0;
	FounderAlleleGraph fag(ped, map);
    
	for(unsigned i = 0; i < map->num_markers(); ++i) {
		fag.reset();
        
        if(not fag.populate(*this, i)) {
            return LOG_ZERO;
        }
		
		if(not fag.likelihood(&tmp_prob, i)) {
			return LOG_ZERO;
        }
		
		return_prob += fag.descentstate_likelihood(i);
    }
	
	return return_prob;
}

string DescentGraph::debug_string() {
    stringstream ss;
	
    ss << "DescentGraph: ";
    
    for(unsigned locus = 0; locus < map->num_markers(); ++locus) {
        for(unsigned i = 0; i  < ped->num_members(); ++i) {
            if(not (ped->get_by_index(i))->isfounder()) {
                ss  << int(get(i, locus, MATERNAL)) \
                    << int(get(i, locus, PATERNAL));
            }
        }
        ss << " ";
    }
    //ss << "\n";
    
    return ss.str();
}
