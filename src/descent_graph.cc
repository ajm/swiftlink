using namespace std;

#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "descent_graph.h"
#include "descent_graph_diff.h"
#include "descent_graph_types.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "genotype.h"
#include "founder_allele_graph.h"
#include "elimination.h"


DescentGraph::DescentGraph(Pedigree* ped, GeneticMap* map) : 
    data(NULL),
    ped(ped), 
    map(map), 
    prob(0.0), 
    marker_transmission(log(0.5) * (2 * (ped->num_members() - ped->num_founders()))), 
    sum_prior_probs(NULL),
    graph_size(2 * ped->num_members()),
    recombinations(-1) {
    
//	graph_size = 2 * ped->num_members();

    //marker_transmission = log(0.5) * 
	//	(ped->num_markers() * 2 * (ped->num_members() - ped->num_founders()));

	// ajm: don't bother calculating the whole likelihood as this is just a 
	// constant for a given pedigree
//	marker_transmission = log(0.5) * 
//		(2 * (ped->num_members() - ped->num_founders()));
	
    data = new char[graph_size * ped->num_markers()];
    
    for(unsigned i = 0; i < graph_size * ped->num_markers(); ++i) {
        data[i] = 0;
    }
    
    sum_prior_probs = new double[ped->num_markers()];
}

DescentGraph::DescentGraph(const DescentGraph& d) : 
    data(NULL),
    ped(d.ped), 
    map(d.map), 
    prob(d.prob), 
	marker_transmission(d.marker_transmission),
    sum_prior_probs(NULL),
    graph_size(d.graph_size),
    recombinations(d.recombinations) {

    unsigned int data_length = graph_size * ped->num_markers();
    
	data = new char[data_length];
    copy(d.data, 
         d.data + data_length, 
         data);
    
    sum_prior_probs = new double[ped->num_markers()];
    copy(d.sum_prior_probs, 
         d.sum_prior_probs + ped->num_markers(), 
         sum_prior_probs);
}

DescentGraph::~DescentGraph() {
	delete[] data;
	delete[] sum_prior_probs;
}

DescentGraph& DescentGraph::operator=(const DescentGraph& d) {
    
	if(&d != this) {
		ped = d.ped;
		map = d.map;
		prob = d.prob;
		graph_size = d.graph_size;
		marker_transmission = d.marker_transmission;
        
        //delete[] data;

        unsigned int data_length = graph_size * ped->num_markers();
    
    	//data = new char[data_length];
        copy(d.data, 
             d.data + data_length, 
             data);
        
        copy(d.sum_prior_probs, \
             d.sum_prior_probs + ped->num_markers(), \
             sum_prior_probs);
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
	
    *p = prob = trans + prior;
    
    //printf("trans = %e, sum prior = %e\n", trans, prior);
	
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
    double theta;
    double antitheta;
	enum parentage parent;
	bool crossover;
	Person* p;
	
    // i = member of family
    // j = mother and father of i
    // k = current loci being considered
    
    recombinations = 0;
	
	for(unsigned k = 0; k < (ped->num_markers() - 1); ++k) { // every loci
	
	    theta = map->get_theta(k);
	    antitheta = map->get_inverse_theta(k);
	
        for(unsigned i = 0; i < ped->num_members(); ++i) { // every person
	    	p = ped->get_by_index(i);

            if(p->isfounder())
                continue;
        
            for(unsigned j = 0; j < 2; ++j) { // mother and father
                parent = static_cast<enum parentage>(j);
                crossover = get(i, k, parent) != get(i, k+1, parent);
            
                if(crossover)
                    recombinations++;
				
                tmp += crossover ? theta : antitheta ;
            }
        }
    }
    
    return tmp;
}

double DescentGraph::_recombination_prob2(unsigned locus) {
    double tmp = 0.0;
    double theta;
    double antitheta;
	enum parentage parent;
	bool crossover;
	Person* p;
	
    // i = member of family
    // j = mother and father of i
    // k = current loci being considered
    
    recombinations = 0;
	
	//for(unsigned k = 0; k < (ped->num_markers() - 1); ++k) { // every loci
	    unsigned k = locus;
	
	    theta = map->get_theta(k);
	    antitheta = map->get_inverse_theta(k);
	    
        for(unsigned i = 0; i < ped->num_members(); ++i) { // every person
	    	p = ped->get_by_index(i);

            if(p->isfounder())
                continue;
        
            for(unsigned j = 0; j < 2; ++j) { // mother and father
                parent = static_cast<enum parentage>(j);
                crossover = get(i, k, parent) != get(i, k+1, parent);
            
                if(crossover)
                    recombinations++;
				
                tmp += crossover ? theta : antitheta ;
            }
        }
    //}
    
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
        
//        printf("%e\n", tmp_prob);
        //fag.print();
        //abort();
        
        sum_prior_probs[i] = tmp_prob;
		
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
		fag.reset();
        
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
	
    //fprintf(stdout, "DescentGraph: \n");
    fprintf(stdout, "DescentGraph: ");
    
    for(int locus = 0; locus < int(ped->num_markers()); ++locus) {
        for(unsigned i = 0; i  < ped->num_members(); ++i) {
            p = ped->get_by_index(i);
            
            if(not p->isfounder()) {
                mat = get(i, locus, MATERNAL);
                pat = get(i, locus, PATERNAL);
                
                //fprintf(stdout, "%d %d %d%d\n", i, locus, mat, pat);
                fprintf(stdout, "%d%d", mat, pat);
            }
        }
        printf(" ");
        //printf("\n");
    }
    
    //fprintf(stdout, "\n");
}

// this seems to be what simwalk2 gets for 'transmission-of-the-markers'
// line ~37200
double DescentGraph::trans_prob() {
    //return marker_transmission;
    return marker_transmission + _recombination_prob();
}

char DescentGraph::get_opposite(unsigned person_id, unsigned locus, enum parentage p) {
    return data[_offset(person_id, locus, p)] == 0 ? 1 : 0;
}

void DescentGraph::write_diff(DescentGraphDiff& dgd) {
    
    for(unsigned i = 0; i < dgd.num_transitions(); ++i) {
        Transition& t = dgd.get_transition(i);
        
        flip_bit(t.get_person(), 
                 t.get_locus(), 
                 t.get_parent());
    }
}

bool DescentGraph::evaluate_diff_sumprior(DescentGraphDiff& diff, double* new_prob) {
    bool ret = true;
    unsigned locus = diff.get_transition(0).get_locus();
    double diff_sumprior_prob;
    
    
    write_diff(diff);
    
    FounderAlleleGraph fag(map, ped);
    fag.reset();
    
    if(not fag.populate(*this, locus)) {
        ret = false;
        goto evaluate_diff_sumprior_end;
    }
    
    if(not fag.likelihood(&diff_sumprior_prob, locus)) {
        ret = false;
        goto evaluate_diff_sumprior_end;
    }
    
    *new_prob = diff_sumprior_prob - sum_prior_probs[locus];
    
 evaluate_diff_sumprior_end:
 
    write_diff(diff);
    
    return ret;
}

void DescentGraph::evaluate_diff_transmission(DescentGraphDiff& diff, double* new_prob, int* new_recombinations) {

    int recombinations_diff = 0;
    double diff_prob = 0.0;
    double orig_prob = 0.0;
    unsigned person;
    unsigned locus;
    enum parentage parent;
    char value;
    
    for(unsigned i = 0; i < diff.num_transitions(); ++i) {
        Transition& t = diff.get_transition(i);
        
        person = t.get_person();
        locus =  t.get_locus();
        parent = t.get_parent();
        
        value = get_opposite(person, locus, parent);
        
        // update transmission probability
        // on the left
        if(locus != 0) {
            if(get(person, locus-1, parent) == value) {
                diff_prob += map->get_inverse_theta(locus-1);
                orig_prob += map->get_theta(locus-1);
                recombinations_diff--;
            }
            else {
                diff_prob += map->get_theta(locus-1);
                orig_prob += map->get_inverse_theta(locus-1);
                recombinations_diff++;
            }
        }
        
        // on the right
        if(locus != (map->num_markers() - 1)) {
            if(get(person, locus+1, parent) == value) {
                diff_prob += map->get_inverse_theta(locus);
                orig_prob += map->get_theta(locus);
                recombinations_diff--;
            }
            else {
                diff_prob += map->get_theta(locus);
                orig_prob += map->get_inverse_theta(locus);
                recombinations_diff++;
            }
        }
    
    }
    
    *new_prob = diff_prob - orig_prob;
    *new_recombinations = recombinations_diff;
}

bool DescentGraph::evaluate_diff(DescentGraphDiff& diff, double* new_prob) {
    
    int recombinations_diff;
    double transmission_diff;
    double sumprior_diff;
    
    evaluate_diff_transmission(diff, &transmission_diff, &recombinations_diff);

    if(not evaluate_diff_sumprior(diff, &sumprior_diff)) {
        return false;
    }
        
    *new_prob = prob + transmission_diff + sumprior_diff;
    
    diff.set_sumprior(sumprior_diff);
    diff.set_prob(*new_prob);
    diff.set_recombinations(recombinations_diff);
    
    return true;
}

void DescentGraph::apply_diff(DescentGraphDiff& diff) {
    
    unsigned locus = diff.get_transition(0).get_locus();
    
    write_diff(diff);
    
    sum_prior_probs[locus] = diff.get_sumprior();
    prob                   = diff.get_prob();
    recombinations += diff.get_recombinations();
}

