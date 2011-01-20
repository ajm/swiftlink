using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

#include "geneticalgorithm.h"
#include "descentgraph.h"
#include "statisticswriter.h"
#include "pedigree.h"

struct compareDescentGraphsPtr {
  bool operator()(const DescentGraph *lhs, const DescentGraph *rhs) { return ! (*lhs < *rhs); }
};

struct compareDescentGraphs {
    bool operator()(const DescentGraph& lhs, const DescentGraph& rhs) { return (lhs.get_prob() > rhs.get_prob()); }
};


// allocs a DescentGraph
void GeneticAlgorithm::crossover_single(DescentGraph* p1, DescentGraph* p2, 
                                        DescentGraph* c1, DescentGraph* c2) {
    unsigned len = p1->data_length();
    unsigned i = rand() % len;

	*c1 = *p1;
    *c2 = *p2;
	
    c1->copy_from(*p2, i, len);
    c2->copy_from(*p1, i, len);
}

void GeneticAlgorithm::crossover_double(DescentGraph* p1, DescentGraph* p2, 
                                        DescentGraph* c1, DescentGraph* c2) {
    unsigned len = p1->data_length();
    unsigned i = rand() % len;
    unsigned j;
    
    while((j = rand() % len) == i)
        ;
    
    unsigned min = i < j ? i : j;
    unsigned max = i < j ? j : i;
    
	*c1 = *p1;
    *c2 = *p2;
    
    c1->copy_from(*p2, min, max);
    c2->copy_from(*p1, min, max);
}

void GeneticAlgorithm::crossover_uniform(DescentGraph* p1, DescentGraph* p2, 
                                         DescentGraph* c1, DescentGraph* c2) {
    unsigned len = p1->data_length();
    double p = 2 / double(len); // on average 2 recombinations
    char tmp1, tmp2;
    
    *c1 = *p1;
    *c2 = *p2;
    
    for(unsigned i = 0; i < len; ++i) {
        if((rand() / double(RAND_MAX)) < p) {
            tmp1 = c1->get_bit(i);
            tmp2 = c2->get_bit(i);
            
            c1->set_bit(i, tmp2);
            c2->set_bit(i, tmp1);
        }
    }
}

// always sorted into descending order
DescentGraph* GeneticAlgorithm::truncation_selection(vector<DescentGraph>& pop) {
    return &pop[rand() % truncation];
}

DescentGraph* GeneticAlgorithm::probabilistic_tournament_selection(
                                    vector<DescentGraph>& pop, unsigned t, double prob) {
    DescentGraph* p[t];
        
	for(unsigned i = 0; i < t; ++i) {
        p[i] = &pop[rand() % pop.size()];
    }

	sort(p, p+t, compareDescentGraphsPtr());
    
	for(unsigned i = 0; i < t; ++i) {
		if((rand() / double(RAND_MAX)) < prob) {
            return p[i];
        }
	}
	
	return p[0];
}

DescentGraph* GeneticAlgorithm::deterministic_tournament_selection(
                                    vector<DescentGraph>& pop, unsigned t) {
	DescentGraph* p[t];
    
	for(unsigned i = 0; i < t; ++i) {
        p[i] = &pop[rand() % pop.size()];
    }
    
	sort(p, p+t, compareDescentGraphsPtr());
    
	return p[0];
}

// this fails completely, the difference between the best and second best
// solutions is seven orders of magnitude
DescentGraph* GeneticAlgorithm::fitness_proportionate_selection(
                                    vector<DescentGraph>& pop, vector<double>& cdf) {
	double i = log(rand() / double(RAND_MAX)) + cdf[cdf.size()-1];

/*
	printf("rand = %f\n", i);
    printf("cdf = %f --> %f\n", cdf[0], cdf[cdf.size()-1]);
    exit(0);
*/

	// quick linear search to test	
    for(unsigned j = 1; j < cdf.size(); ++j)
        if((cdf[j-1] < i) and (i <= cdf[j]))
            return &pop[j];
    
    return &pop[pop.size()-1];
}

// mutate c, and return it
DescentGraph* GeneticAlgorithm::mutate(DescentGraph* c, double prob) {
    unsigned len = c->data_length();
    //double p = 1 / double(len);
    
    for(unsigned i = 0; i < len; ++i) {
        if((rand() / double(RAND_MAX)) < prob) {
            c->flip_bit(i);
        }
    }
    
    return c;
}

unsigned GeneticAlgorithm::get_random(unsigned i) {
	return rand() % i;
}

unsigned GeneticAlgorithm::get_random_nonfounder() {
	return get_random(ped->num_members() - ped->num_founders()) + ped->num_founders();
}

unsigned GeneticAlgorithm::get_random_locus() {
	return get_random(ped->num_markers());
}

enum parentage GeneticAlgorithm::get_random_parent() {
	return static_cast<enum parentage>(get_random(2));
}

// mutate c, and return it
DescentGraph* GeneticAlgorithm::mutate(DescentGraph* c, int freq) {
    
    for(int i = 0; i < freq; ++i) {
        c->flip_bit(get_random_nonfounder(), \
                    get_random_locus(), \
                    get_random_parent());
    }
    
    return c;
}


void GeneticAlgorithm::make_cdf(vector<DescentGraph>& pop, vector<double>& cdf) {
	cdf[0] = pop[0].get_prob();
        
	for(unsigned i = 1; i < pop.size(); ++i) {
			
		if(pop[i].illegal()) {
			cdf[i] = cdf[i-1];
			continue;
		}
            
		double tmplik = pop[i].get_prob();
			
		cdf[i] = log(exp(cdf[i-1] - tmplik) + 1) + tmplik;
	}
}

// TODO: figure out when termination happens early
//DescentGraph* GeneticAlgorithm::optimise(unsigned popsize, unsigned iterations, 
//                                         unsigned elitism) {
DescentGraph* GeneticAlgorithm::optimise(unsigned popsize, unsigned iterations, unsigned elitism, 
		                  int crossover, int selection, unsigned selection_size, double selection_prob) {
    DescentGraph *pa, *pb, *tmp;
	DescentGraph *best;
    double mutation_prob;
    
	vector<DescentGraph> pop(popsize, DescentGraph(ped, map));
    vector<DescentGraph> newpop(popsize, DescentGraph(ped, map));
    //vector<double> cdf(popsize, 0.0);
    
    StatisticsWriter sw;
    double values[popsize];
    
    
	best = new DescentGraph(ped, map);
	best->random();
	//best->likelihood();
	best->haplotype_likelihood();
    
    mutation_prob = 1 / double(best->data_length());
    truncation = popsize * selection_prob;
    
    fprintf(stderr, "population = %u\n", popsize);
    fprintf(stderr, "iterations = %u\n", iterations);
    fprintf(stderr, "elitism = %u\n", elitism);
    fprintf(stderr, "crossover = %d\n", crossover);
    fprintf(stderr, "selection = %d\n", selection);
    fprintf(stderr, "selection size = %u\n", selection_size);
    fprintf(stderr, "selection prob = %f\n", selection_prob);
    
    
    // initialise each element
    for(unsigned i = 0; i < pop.size(); ++i) {
        // genotype elimination should be done separately
        // and then passed to DescentGraph
        if(not pop[i].random()) {
            abort();
        }
    }

	fprintf(stderr, "finished init\n");
    
    // run genetic algorithm
    for(unsigned i = 0; i < iterations; ++i) {
        // assess fitness
        for(unsigned j = 0; j < pop.size(); ++j) {
            //pop[j].likelihood();
			//pop[j].likelihood(&values[j]);
            //pop[j].haplotype_likelihood();
            pop[j].haplotype_likelihood(&values[j]);
//            pop[j].print();
        }
        
//        return best;
        
        sw.print(values, popsize);
                
		// needed for truncation selection & best tracking
        sort(pop.begin(), pop.end(), compareDescentGraphs());
        
		// needed for fitness proportional selection
        //make_cdf(pop, cdf);

/*
		printf("Population:\n");
		for(unsigned j = 0; j < pop.size(); ++j) {
            //printf("\t%f\n", pop[j].get_prob() / log(10));
            pop[j].print();
        }
*/

		// is there a new best?
		tmp = &pop[0];
		if(tmp->get_prob() > best->get_prob()) {
			*best = *tmp;
		}

		if((i % 10) == 0)
			printf("%u %f\n", i, best->get_prob() / log(10));
        
        // put fittest solutions from the current population
        // into the next population
        for(unsigned j = 0; j < elitism; ++j) {
            newpop[j] = pop[j];
        }
        
        // create new population
   	    for(unsigned j = elitism / 2; j < (pop.size() / 2); ++j) {
   	        //pa = truncation_selection(pop);
   	        //pb = truncation_selection(pop);
            
			//pa = fitness_proportionate_selection(pop, cdf);
			//pb = fitness_proportionate_selection(pop, cdf);

            //pa = deterministic_tournament_selection(pop, 2);
   	        //pb = deterministic_tournament_selection(pop, 2);
            
            //pa = probabilistic_tournament_selection(pop, 5, 0.25);
   	        //pb = probabilistic_tournament_selection(pop, 5, 0.25);
   	        
   	        switch(selection) {
   	            case 0 :
   	                pa = truncation_selection(pop);
   	                pb = truncation_selection(pop);
   	                break;
   	                
   	            case 1 :
   	                pa = probabilistic_tournament_selection(pop, selection_size, selection_prob);
   	                pb = probabilistic_tournament_selection(pop, selection_size, selection_prob);
   	                break;
   	            
   	            default :
   	                abort();
   	        }
   	        
   	        if((rand() / double(RAND_MAX)) < 0.5) {
            
            switch(crossover) {
                case 0:
                    crossover_single(pa, pb, &newpop[2*j], &newpop[(2*j)+1]);
                    break;

                case 1:
                    crossover_double(pa, pb, &newpop[2*j], &newpop[(2*j)+1]);
                    break;

                case 2:
                    crossover_uniform(pa, pb, &newpop[2*j], &newpop[(2*j)+1]);
                    break;
                
                default :
                    abort();
            }
            }
            else {
/*
   	        crossover_single(
			//crossover_double(
			//crossover_uniform(
                             pa, 
                             pb, 
                             &newpop[2*j], 
                             &newpop[(2*j)+1] );
*/   	    
            //mutate(&newpop[2*j], mutation_prob);
            //mutate(&newpop[(2*j)+1], mutation_prob);
            
            newpop[2*j] = *pa;
            newpop[(2*j)+1] = *pb;
            
            mutate(&newpop[2*j], 1);
            mutate(&newpop[(2*j)+1], 1);
            }
        }
        
		pop.swap(newpop);
    }
    
    return best;
}

