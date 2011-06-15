using namespace std;

#include <cstdio>
#include <string>
#include <fstream>

#include "haplotype_writer.h"
#include "descent_graph.h"
#include "founder_allele_graph.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "person.h"


// this is how genehunter writes out haplotypes
bool HaplotypeWriter::write() {
	fstream hap;
	Person* p;
	//int paternal_alleles[ped->num_members()][ped->num_markers()]; // :-(
    //int maternal_alleles[ped->num_members()][ped->num_markers()];
    int* paternal_alleles;
    int* maternal_alleles;
    FounderAlleleGraph fag(map, ped);
    double prob;
    
    paternal_alleles = new int[ped->num_members() * ped->num_markers()];
    maternal_alleles = new int[ped->num_members() * ped->num_markers()];
    
    // fill in alleles
    for(unsigned j = 0; j < ped->num_markers(); ++j) {
        fag.reset();
        fag.populate(*dg, j);
        fag.likelihood(&prob, j);
        
        for(unsigned i = 0; i < ped->num_members(); ++i) {
            //paternal_alleles[i][j] = fag.get_founderallele_assignment(dg->get_founderallele(i, j, MATERNAL));
            //maternal_alleles[i][j] = fag.get_founderallele_assignment(dg->get_founderallele(i, j, PATERNAL));
            paternal_alleles[(ped->num_members() * i) + j] = fag.get_founderallele_assignment(dg->get_founderallele(i, j, MATERNAL));
            maternal_alleles[(ped->num_members() * i) + j] = fag.get_founderallele_assignment(dg->get_founderallele(i, j, PATERNAL));
        }
    }

    fprintf(stderr, 
            "Printing haplotype (log-likelihood %.5f) to %s\n", 
            dg->get_prob() / log(10), filename.c_str());
    
    // actually print to file
	hap.open(filename.c_str(), ios::out);
    
	if(not hap.is_open()) {
		fprintf(stderr, "Could not open haplotype output file \"%s\"\n", 
			filename.c_str());
		return false;
	}
    
	for(unsigned i = 0; i < ped->num_members(); ++i) {
		p = ped->get_by_index(i);

		hap << p->get_id() << '\t' \
			<< p->get_father() << '\t' \
			<< p->get_mother() << '\t' \
			<< int(p->get_affection()) << '\t' ;

		for(unsigned j = 0; j < p->num_markers(); ++j) {
			hap << ' ' << paternal_alleles[(ped->num_members() * i) + j];
		}
		hap << endl;


		hap << "\t\t\t\t";

		for(unsigned j = 0; j < p->num_markers(); ++j) {
			hap << ' ' << maternal_alleles[(ped->num_members() * i) + j];
		}
		hap << endl;
	}
    
	hap.close();
    
    delete[] paternal_alleles;
    delete[] maternal_alleles;
	
	return true;
}

