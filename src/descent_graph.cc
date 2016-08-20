#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

#include "types.h"
#include "descent_graph.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "elimination.h"
#include "founder_allele_graph4.h"

using namespace std;


DescentGraph::DescentGraph(Pedigree* ped, GeneticMap* map, bool sex_linked) :
    data(NULL),
    ped(ped), 
    map(map), 
    prob(0.0), 
    marker_transmission(log(0.5) * (2 * (ped->num_members() - ped->num_founders()))), 
    graph_size(2 * ped->num_members()),
    recombinations(-1),
    sex_linked(sex_linked),
    seq() {
    
    data = new int[graph_size * map->num_markers()];
    
    for(unsigned i = 0; i < (graph_size * map->num_markers()); ++i) {
        data[i] = 0;
    }
    
    if(sex_linked) {
        marker_transmission = log(0.5) * (ped->num_members() - ped->num_founders());
        //_invalidate_paternal_x();
    }

    find_founderallelegraph_ordering();
}

DescentGraph::DescentGraph(const DescentGraph& d) : 
    data(NULL),
    ped(d.ped), 
    map(d.map), 
    prob(d.prob), 
	marker_transmission(d.marker_transmission),
    graph_size(d.graph_size),
    recombinations(d.recombinations),
    sex_linked(d.sex_linked),
    seq(d.seq) {

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
        sex_linked = d.sex_linked;
                
        copy(d.data, 
             d.data + (graph_size * map->num_markers()), 
             data);
    }

	return *this;
}

void DescentGraph::_invalidate_paternal_x() {
    Person* p;

    for(unsigned i = 0; i < ped->num_members(); ++i) { // every person
        p = ped->get_by_index(i);
        
        if(p->isfounder())
            continue;

        if(p->ismale()) {
            for(unsigned j = 0; j < map->num_markers(); ++j) { // every loci
                set(i, j, PATERNAL, NONE);
            }
        }
    }
}

bool DescentGraph::random_descentgraph() {
    GenotypeElimination ge(ped, sex_linked);
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

void DescentGraph::copy_locus(unsigned src, unsigned dst) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        data[_offset(i, dst, MATERNAL)] = data[_offset(i, src, MATERNAL)];
        data[_offset(i, dst, PATERNAL)] = data[_offset(i, src, PATERNAL)];
    }
}

void DescentGraph::copy_locus(DescentGraph& d, unsigned src, unsigned dst) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        data[_offset(i, dst, MATERNAL)] = d.get(i, src, MATERNAL);
        data[_offset(i, dst, PATERNAL)] = d.get(i, src, PATERNAL);
    }
}

void DescentGraph::set(unsigned person_id, unsigned locus, enum parentage p, int value) {
    data[_offset(person_id, locus, p)] = value;
}

void DescentGraph::flip(unsigned person_id, unsigned locus, enum parentage p) {
	flip_bit(_offset(person_id, locus, p));
}

// this only works because the founders are guaranteed to be at the start of
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
    //fprintf(stderr, "TP = %f, SPP = %f\n", _transmission_prob(), _sum_prior_prob());
    //prob = _transmission_prob();
    return prob;
}
/*
double DescentGraph::get_haplotype_likelihood() {
    prob = log_product(_transmission_prob(), _best_prior_prob());
    return prob;
}
*/
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

    unsigned num_alleles = sex_linked ? 1 : 2;

    for(unsigned i = 0; i < ped->num_members(); ++i) { // every person
        p = ped->get_by_index(i);

        if(p->isfounder())
            continue;

        for(unsigned j = 0; j < num_alleles; ++j) { // mother and father, just mother (0) for sex-linked
            parent = static_cast<enum parentage>(j);
            crossover = get(i, locus, parent) != get(i, locus + 1, parent);

            if(count_crossovers and crossover) {
                ++recombinations;
            }

            tmp += (crossover ? theta : antitheta) ;
        }
    }
    
    return tmp;
}

double DescentGraph::_sum_prior_prob() {
    double tmp_prob;
	double return_prob = 0.0;
    FounderAlleleGraph4 f(ped, map, sex_linked);
    
    f.set_sequence(&seq);
    
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        f.set_locus(i);
        f.reset(*this);
        
        if((tmp_prob = f.likelihood()) == 0.0) {
            fprintf(stderr, "error: descent graph illegal at locus %d\n", int(i));
            fprintf(stderr, "%s\n", f.debug_string().c_str());
            fprintf(stderr, "%s\n", debug_string().c_str());
			return LOG_ZERO;
        }
        
		return_prob += log(tmp_prob);
    }
    
    return return_prob;
    
}

/*
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
*/

string DescentGraph::debug_string() {
    stringstream ss;
    Person* p;
	
    ss << "DescentGraph order: ";
    
    for(unsigned i = 0; i  < ped->num_members(); ++i) {
        p = ped->get_by_index(i);
        if(not p->isfounder()) {
            ss << p->get_id() << " ";
        }
    }
    ss << "\n";

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

void DescentGraph::find_founderallelegraph_ordering() {
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

