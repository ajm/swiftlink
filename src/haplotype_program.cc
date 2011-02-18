using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "linkageprogram.h"
#include "linkagefileparser.h"
#include "pedfileparser.h"
#include "mapfileparser.h"
#include "geneticmap.h"
#include "diseasemodel.h"
#include "descentgraph.h"
#include "geneticalgorithm.h"
#include "simulatedannealing.h"
#include "haplotypewriter.h"
#include "configparser.h"


bool HaplotypeProgram::run() {
    return run_ga();
}

bool HaplotypeProgram::run_ga() {
	DescentGraph* ans;

	if(! read_and_check_input()) {
		return false;
	}
/*    
    for(int i = 0; i < int(pedigrees.size()); ++i)
		pedigrees[i]->print();
	
	map.print();
	dm.print();
*/    
    srand(time(NULL));
    
    ConfigParser cp("config.txt");
    if(! cp.parse())
		return false;
//    unsigned population = 100;
//    unsigned iterations = 1000;
//    unsigned elitism = 0;
	
	for(int i = 0; i < int(pedigrees.size()); ++i) {
		GeneticAlgorithm ga(pedigrees[i], &map);
		ans = ga.optimise(cp.population, cp.iterations, cp.elitism, 
		                  cp.crossover, cp.selection, cp.selection_size, cp.selection_prob);
		
		//printf("result = %f\n", ans->get_prob() / log(10));
	
		HaplotypeWriter hw("haplotype.txt", ans, pedigrees[i], &map);
		hw.write();
		
		delete ans;
	}
	
	return true;
}

bool HaplotypeProgram::run_sa() {
	SA_DescentGraph* ans;

	if(! read_and_check_input()) {
		return false;
	}

	srand(time(NULL));

	ConfigParser cp("config.txt");
    if(! cp.parse())
		return false;

	for(int i = 0; i < int(pedigrees.size()); ++i) {
		SimulatedAnnealing sa(pedigrees[i], &map);
		ans = sa.optimise(1600 * pedigrees[i]->num_members() * pedigrees[i]->num_markers() * 20 * 2);

		HaplotypeWriter hw("haplotype.txt", ans, pedigrees[i], &map);
		hw.write();
		
		delete ans;
	}

	return true;
}

