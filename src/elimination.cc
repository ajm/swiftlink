using namespace std;

#include <cstdlib>

#include "elimination.h"
#include "pedigree.h"
#include "person.h"
#include "descent_graph.h"
#include "genotype.h"

void GenotypeElimination::_initial_elimination() {
    Person* tmp;
    enum unphased_genotype g;
    
	for(unsigned i = 0; i < ped->num_markers(); ++i) {		
		// all possible initially at this locus
		for(unsigned j = 0; j < ped->num_members(); ++j) {
			genotype_set_all(possible_genotypes[i][j]);
        }
		
		// list compatible genotypes
		// (ie: mark which phased genotypes are impossible based on unphased genotypes)
		for(unsigned j = 0; j < ped->num_members(); ++j) {
			tmp = ped->get_by_index(j);
			
			g = tmp->get_genotype(i);
			
			switch(g) {
				case HOMOZ_A:
					genotype_set_homozA(possible_genotypes[i][j]);
					break;
					
				case HETERO:
					genotype_set_hetero(possible_genotypes[i][j]);
					break;
					
				case HOMOZ_B:
					genotype_set_homozB(possible_genotypes[i][j]);
					break;
					
				case UNTYPED:
					genotype_set_all(possible_genotypes[i][j]);
					break;
					
				default :
					break;
			}
		}
    }
}

bool GenotypeElimination::_elimination_pass(int** ds, unsigned locus) {
    Person* tmp;
    int changes, mat, pat;
    
    while(true) {
        changes = 0;
			
        for(unsigned i = 0; i < ped->num_members(); ++i) {
            tmp = ped->get_by_index(i);
				
            if(tmp->isfounder())
                continue;
                
            pat = tmp->get_paternalid();
            mat = tmp->get_maternalid();
                
            changes += _parent_homoz(ds, locus, mat, pat, i, AA, MATERNAL);
            changes += _parent_homoz(ds, locus, mat, pat, i, BB, MATERNAL);
            changes += _parent_homoz(ds, locus, pat, mat, i, AA, PATERNAL);
            changes += _parent_homoz(ds, locus, pat, mat, i, BB, PATERNAL);
            
            changes += _child_homoz(ds, locus, mat, pat, i, AA);
            changes += _child_homoz(ds, locus, mat, pat, i, BB);                
        }
        
        if(not _legal(ds, locus)) {
            return false;
        }
			
        if(changes == 0) {
            break;
        }
    }
    
    return true;
}

bool GenotypeElimination::_legal(int** ds, unsigned locus) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        if(genotype_bad(ds[locus][i])) {
			return false;
        }
	}
	
    return true;
}

bool GenotypeElimination::_complete(int** ds, unsigned locus) {
	for(unsigned i = 0; i < ped->num_members(); ++i) {
		if(not genotype_single(ds[locus][i])) {
			return false;
        }
	}
	
	return true;
}

int GenotypeElimination::_child_homoz(int** ds, 
                                      unsigned locus, 
                                      int mother, 
                                      int father, 
                                      int child,
                                      enum phased_genotype homoz) {
	enum phased_genotype other_homoz;
	int changes = 0;
	
	other_homoz = (homoz == AA) ? BB : AA;
	
	if( genotype_only(ds[locus][child], homoz) ) {
		
		if( genotype_possible(ds[locus][mother], other_homoz) ) {
			genotype_remove(ds[locus][mother], other_homoz);
			changes++;
		}
		
		if( genotype_possible(ds[locus][father], other_homoz) ) {
			genotype_remove(ds[locus][father], other_homoz);
			changes++;
		}
	}
	
	return changes;
}

int GenotypeElimination::_parent_homoz(int** ds,
                                       unsigned locus,
                                       int parent, 
                                       int other_parent, 
                                       int child, 
                                       enum phased_genotype homoz, 
                                       enum parentage p) {
	int changes = 0;
	enum phased_genotype other_homoz, 
                         hetero, 
                         other_hetero;
	
	// only homoz parents mean children cannot have other allele
	// from that parent
	if(genotype_only(ds[locus][parent], homoz)) {
		
		if(((p == MATERNAL) && (homoz == AA)) || ((p == PATERNAL) && (homoz == BB))) {
			hetero = AB;
			other_hetero = BA;
		}
		else {
			hetero = BA;
			other_hetero = AB;
		}
		
		other_homoz = (homoz == AA) ? BB : AA;
        
		
		if(genotype_possible(ds[locus][child], other_homoz)) {
			genotype_remove(ds[locus][child], other_homoz);
			changes++;
		}
		
		if(genotype_possible(ds[locus][child], other_hetero)) {
			genotype_remove(ds[locus][child], other_hetero);
			changes++;
		}
		
		// if child hetero only, then other_parent cannot be homoz
		if(genotype_possible(ds[locus][child], hetero) && \
          !genotype_possible(ds[locus][child], homoz) ) {
			
			if(genotype_possible(ds[locus][other_parent], homoz) ) {
                genotype_remove(ds[locus][other_parent], homoz);
				changes++;
			}
		}
	}
	
	return changes;
}

void GenotypeElimination::_random_eliminate(int** ds, unsigned locus) {
	//int queue[ped->num_members()];
    vector<int> queue(ped->num_members());
    int qindex, tmp, tmp2;
    
	// find everyone for whom their genotype is ambiguous
	qindex = 0;
	for(unsigned i = 0; i < ped->num_members(); ++i) {
		if(! genotype_single(ds[locus][i])) {
			queue[qindex++] = i;
		}
	}
	
	// get random person
	tmp = queue[random() % qindex];
	
	// find all possible genotypes for tmp
	qindex = 0;
	for(int i = 0; i < 4; ++i ) {
		tmp2 = 1 << i;
		if(genotype_possible(ds[locus][tmp], tmp2)) {
			queue[qindex++] = tmp2;
        }
	}
	
	tmp2 = queue[random() % qindex];
	
	ds[locus][tmp] = tmp2;
}

enum parentage GenotypeElimination::_state(enum phased_genotype child, 
                      enum phased_genotype parent, 
                      enum parentage p) {
	
    enum phased_genotype hetero;
    
    hetero = (p == MATERNAL) ? AB : BA;
	
	if(child == AA || child == hetero) {
		if(parent == AB) {
			return MATERNAL;
        }
		else if(parent == BA) {
			return PATERNAL;
        }
		else {
			return static_cast<enum parentage>(rand() % 2);
        }
    }
	else {
		if(parent == BA) {
			return MATERNAL;
        }
        else if(parent == AB) {
			return PATERNAL;
        }
		else {
			return static_cast<enum parentage>(rand() % 2);
        }
	}
	
	abort();
}

void GenotypeElimination::_write_descentgraph(DescentGraph& d, int** ds) {
    Person* tmp;
    enum phased_genotype m, p, c;
    
    for(unsigned i = 0; i < ped->num_markers(); ++i) {
        for(unsigned j = 0; j < ped->num_members(); ++j) {
            tmp = ped->get_by_index(j);
            
            if(tmp->isfounder()) {
                d.set(j, i, MATERNAL, 0);
                d.set(j, i, PATERNAL, 0);
                continue;
            }
            
            m = static_cast<enum phased_genotype>(ds[i][tmp->get_maternalid()]);
            p = static_cast<enum phased_genotype>(ds[i][tmp->get_paternalid()]);
            c = static_cast<enum phased_genotype>(ds[i][j]);
            
            d.set(j, i, MATERNAL, _state(c, m, MATERNAL));
            d.set(j, i, PATERNAL, _state(c, p, PATERNAL));
        }
    }
}

void GenotypeElimination::_copy_descentstate(int** src, int** dst) {
    for(unsigned i = 0; i < ped->num_markers(); ++i) {
        _copy_descentstate_locus(src, dst, i);
//        for(unsigned j = 0; j < ped->num_members(); ++j) {
//            dst[i][j] = src[i][j];
//        }
    }
}

void GenotypeElimination::_copy_descentstate_locus(int** src, int** dst, unsigned locus) {
    for(unsigned j = 0; j < ped->num_members(); ++j) {
        dst[locus][j] = src[locus][j];
    }
}

bool GenotypeElimination::random_descentgraph(DescentGraph& d) {
    //int** descentstate = new int**[ped->num_markers()][ped->num_members()];
    int** descentstate = new int*[ped->num_markers()];
    for(int i = 0; i < int(ped->num_markers()); ++i) {
        descentstate[i] = new int[ped->num_members()];
    }
    
    // make sure this is only done once
    if(not init_processing) {
        _initial_elimination();
        for(unsigned i = 0; i < ped->num_markers(); ++i) {
            if(not _elimination_pass(possible_genotypes, i)) {
                return false;
            }
        }
        init_processing = true;
    }
    
    _copy_descentstate(possible_genotypes, descentstate);
    
    // random elimination stuff
    for(unsigned i = 0; i < ped->num_markers(); ++i) {
        while(true) {
            if(not _elimination_pass(descentstate, i)) {
                _copy_descentstate_locus(possible_genotypes, descentstate, i);
            }
            
            if(_complete(descentstate, i)) {
                break;
            }
            else {
                _random_eliminate(descentstate, i);
            }
        }
    }
    
    _write_descentgraph(d, descentstate);
    
    for(int i = 0; i < int(ped->num_markers()); ++i) {
        delete[] descentstate[i];
    }
    delete[] descentstate;
    
    return true;
}
