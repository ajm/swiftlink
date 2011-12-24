using namespace std;

#include <cstdio>
#include <cmath>
#include <queue>
#include <vector>
#include <numeric>
#include <algorithm>

#include "types.h"
#include "person.h"
#include "pedigree.h"
#include "founder_allele_graph4.h"
#include "descent_graph.h"
#include "genetic_map.h"


string FounderAlleleGraph4::debug_string() { 
    stringstream ss;
        
    // ???
    
    return ss.str();
}


// -----------------------------------------
double FounderAlleleGraph4::init_likelihood(DescentGraph& dg, int newlocus) {
    Person* p;
    int pid;
    int parent_allele;
    int tmp;
    
    locus = newlocus;
    major_freq = map->get_major(locus);
    minor_freq = map->get_minor(locus);
    
	// find founder allele assignments, this is only related to the current 
	// descent graph and not whether people are typed or not
	for(unsigned i = 0; i < ped->num_members(); ++i) {
	    pid = (*sequence)[i];
        p = ped->get_by_index(pid);
	    
	    tmp = pid * 2;
	    
	    if(p->isfounder()) {
	        edge_list[tmp] = tmp;
	        edge_list[tmp + 1] = tmp + 1;
	    }
	    else {
	        parent_allele = dg.get(pid, locus, MATERNAL);
	        edge_list[tmp] = edge_list[(p->get_maternalid() * 2) + parent_allele];
	        
	        parent_allele = dg.get(pid, locus, PATERNAL);
	        edge_list[tmp + 1] = edge_list[(p->get_paternalid() * 2) + parent_allele];
	    }
	}
	
	
	enum unphased_genotype g;
    int mat_fa;
    int pat_fa;
    int num_groups;
    int group_index;
    int group1, group2;
    int fixed1, fixed2;
    bool legal0, legal1, legal2, legal3;
    
    // reset everything
    num_groups = 0;
    group_index = 0;
    group_membership.assign(founder_alleles, DEFAULT_COMPONENT);
    group_fixed.assign(founder_alleles, -1);
    allele_assignment[0].assign(founder_alleles, UNTYPED);
    allele_assignment[1].assign(founder_alleles, UNTYPED);
    
    // calculate likelihood, the founder allele graph is only
    // constrained by typed members of the pedigree
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        p = ped->get_by_index(i);
        if(not p->istyped())
            continue;
        
        g = p->get_marker(locus);
        
        if(g == UNTYPED)
            continue;
        
        tmp = i * 2;
        mat_fa = edge_list[tmp];
        pat_fa = edge_list[++tmp];
        
        
        // autozygous
        if(mat_fa == pat_fa) { /* 301 - 347 */
            if(g == HETERO) {
                //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                return 0.0;
            }
            
            group1 = group_membership[mat_fa];
            
            if(group1 != DEFAULT_COMPONENT) {
                // already belongs to a component
                fixed1 = group_fixed[group1];
                if(fixed1 != -1) {
                    if(g != allele_assignment[fixed1][mat_fa]) {
                        //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                        return 0.0;
                    }
                }
                else {
                    if(allele_assignment[0][mat_fa] == g) {
                        group_fixed[group1] = 0;
                    }
                    else if(allele_assignment[1][mat_fa] == g) {
                        group_fixed[group1] = 1;
                    }
                    else {
                        //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                        return 0.0;
                    }
                }
            }
            else {
                group_membership[mat_fa] = group_index;
                group_fixed[group_index] = 0;
                allele_assignment[0][mat_fa] = g;
                ++group_index;
                ++num_groups;
            }
        }
        // not autozygous
        else {
            group1 = group_membership[mat_fa];
            group2 = group_membership[pat_fa];
            
            if(group1 != DEFAULT_COMPONENT) { /* 355 - 548 */
                
                fixed1 = group_fixed[group1];
                
                if(group2 != DEFAULT_COMPONENT) {
                    // both in the same group
                    if(group1 == group2) {
                        if(fixed1 != -1) {
                            // fixed, check if legit
                            if(not legal(g, allele_assignment[fixed1][mat_fa], allele_assignment[fixed1][pat_fa])) {
                                //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                return 0.0;
                            }
                        }
                        else {
                            // not fixed, check if still the case
                            legal0 = legal(g, allele_assignment[0][mat_fa], allele_assignment[0][pat_fa]);
                            legal1 = legal(g, allele_assignment[1][mat_fa], allele_assignment[1][pat_fa]);
                            
                            if(legal0 and legal1) {
                                group_fixed[group1] = -1;
                            }
                            else {
                                if(legal0) {
                                    group_fixed[group1] = 0;
                                }
                                else if(legal1) {
                                    group_fixed[group1] = 1;
                                }
                                else {
                                    //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                    return 0.0;
                                }
                            }
                        }
                    }
                    else { /* 388 - 502 */
                        // in different groups
                        fixed2 = group_fixed[group2];
                        
                        if(fixed1 != -1) {
                            if(fixed2 != -1) {
                                // fixed, check if legit
                                if(not legal(g, allele_assignment[fixed1][mat_fa], allele_assignment[fixed2][pat_fa])) {
                                    //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                    return 0.0;
                                }
                            }
                            else {
                                // group1 is fixed, which assignment for group 2 works
                                legal0 = legal(g, allele_assignment[fixed1][mat_fa], allele_assignment[0][pat_fa]);
                                legal1 = legal(g, allele_assignment[fixed1][mat_fa], allele_assignment[1][pat_fa]);
                                
                                if(legal0)
                                    fixed2 = 0;
                                else if(legal1)
                                    fixed2 = 1;
                                else {
                                    //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                    return 0.0;
                                }
                            }
                        }
                        else if(fixed2 != -1) {
                            // group2 is fixed, which assignment for group 1 works
                            legal0 = legal(g, allele_assignment[0][mat_fa], allele_assignment[fixed2][pat_fa]);
                            legal1 = legal(g, allele_assignment[1][mat_fa], allele_assignment[fixed2][pat_fa]);
                            
                            if(legal0)
                                fixed1 = 0;
                            else if(legal1)
                                fixed1 = 1;
                            else {
                                //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                return 0.0;
                            }
                        }
                        else {
                            // neither group1 nor group2 are fixed
                            legal0 = legal(g, allele_assignment[0][mat_fa], allele_assignment[0][pat_fa]);
                            legal1 = legal(g, allele_assignment[1][mat_fa], allele_assignment[0][pat_fa]);
                            legal2 = legal(g, allele_assignment[0][mat_fa], allele_assignment[1][pat_fa]);
                            legal3 = legal(g, allele_assignment[1][mat_fa], allele_assignment[1][pat_fa]);
                            
                            if(not (legal0 or legal1 or legal2 or legal3)) {
                                //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                return 0.0;
                            }
                            
                            if(legal0 and not(legal1 or legal2 or legal3)) {
                                // fixed, allele 0, no swap
                                fixed1 = fixed2 = 0;
                            }
                            else if(legal1 and not(legal0 or legal2 or legal3)) {
                                // fixed, swap assignment in one group
                                fixed1 = 1;
                                fixed2 = 0;
                            }
                            else if(legal2 and not(legal0 or legal1 or legal3)) {
                                // fixed, swap assignment in one group
                                fixed1 = 0;
                                fixed2 = 1;
                            }
                            else if(legal3 and not(legal0 or legal1 or legal2)) {
                                // fixed, allele 1, no swap
                                fixed1 = fixed2 = 1;
                            }
                            else if(legal0 and legal3 and not (legal1 or legal2)){
                                // still unfixed, assignments don't need swapping
                                fixed1 = fixed2 = -1;
                            }
                            else if(legal1 and legal2 and not (legal0 or legal3)){
                                // still unfixed, swap assignment in one group
                                fixed1 = fixed2 = -2;
                            }
                        }
                        
                        if(fixed1 != fixed2) {
                            // XXX for now, always keep group1
                            group_fixed[group1] = fixed1;
                            combine_components(group1, group2, true);
                        }
                        else {
                            if(fixed1 == -2) {
                                group_fixed[group1] = -1;
                                combine_components(group1, group2, true);
                            }
                            else {
                                group_fixed[group1] = fixed1;
                                combine_components(group1, group2, false);
                            }
                        }
                        
                        --num_groups;
                    }
                }
                else {
                    if(fixed1 != -1) {
                        enum unphased_genotype tmp = get_other_allele(g, allele_assignment[fixed1][mat_fa]);
                        if(tmp == UNTYPED) {
                            //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                            return 0.0;
                        }
                        else 
                            allele_assignment[fixed1][pat_fa] = tmp;
                    }
                    else {
                        enum unphased_genotype tmp0 = get_other_allele(g, allele_assignment[0][mat_fa]);
                        enum unphased_genotype tmp1 = get_other_allele(g, allele_assignment[1][mat_fa]);
                        
                        if(tmp0 != UNTYPED) {
                            if(tmp1 != UNTYPED) {
                                allele_assignment[0][pat_fa] = tmp0;
                                allele_assignment[1][pat_fa] = tmp1;
                            }
                            else {
                                allele_assignment[0][pat_fa] = tmp0;
                                group_fixed[group1] = 0;
                            }
                        }
                        else {
                            if(tmp1 != UNTYPED) {
                                allele_assignment[1][pat_fa] = tmp1;
                                group_fixed[group1] = 1;
                            }
                            else {
                                //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                                return 0.0;
                            }
                        }
                        
                    }
                    
                    group_membership[pat_fa] = group1;
                }
            }
            else if(group2 != DEFAULT_COMPONENT) { /* 550 - 594 */
                fixed2 = group_fixed[group2];
                
                if(fixed2 != -1) {
                    enum unphased_genotype tmp = get_other_allele(g, allele_assignment[fixed2][pat_fa]);
                    if(tmp == UNTYPED) {
                        //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                        return 0.0;
                    }
                    else 
                        allele_assignment[fixed2][mat_fa] = tmp;
                }
                else {
                    enum unphased_genotype tmp0 = get_other_allele(g, allele_assignment[0][pat_fa]);
                    enum unphased_genotype tmp1 = get_other_allele(g, allele_assignment[1][pat_fa]);
                    
                    if(tmp0 != UNTYPED) {
                        if(tmp1 != UNTYPED) {
                            allele_assignment[0][mat_fa] = tmp0;
                            allele_assignment[1][mat_fa] = tmp1;
                        }
                        else {
                            allele_assignment[0][mat_fa] = tmp0;
                            group_fixed[group2] = 0;
                        }
                    }
                    else {
                        if(tmp1 != UNTYPED) {
                            allele_assignment[1][mat_fa] = tmp1;
                            group_fixed[group2] = 1;
                        }
                        else {
                            //fprintf(stderr, "FAG ILLEGAL: %d\n", __LINE__);
                            return 0.0;
                        }
                    }        
                }
                
                group_membership[mat_fa] = group2;            
            }
            else { /* 596 - 621 */
                if(g == HETERO) {
                    allele_assignment[0][mat_fa] = allele_assignment[1][pat_fa] = HOMOZ_A;
                    allele_assignment[1][mat_fa] = allele_assignment[0][pat_fa] = HOMOZ_B;
                    group_fixed[group_index] = -1;
                }
                else {
                    allele_assignment[0][mat_fa] = g;
                    allele_assignment[1][pat_fa] = g;
                    group_fixed[group_index] = 0;
                }
                
                // neither in a group
                group_membership[mat_fa] = group_index;
                group_membership[pat_fa] = group_index;
                ++group_index;
                ++num_groups;
            }
        }
    }
    
    // calculate likelihood
    double prob = 1.0; /* 626 - 642 */
    
    // XXX TODO
    
    return prob;
}

// should be cached, this is temporary
double FounderAlleleGraph4::get_freq(enum unphased_genotype g) {
    return g == HOMOZ_A ? major_freq : minor_freq;
}

// obs is observed genotype (edge), can be anything apart from UNTYPED
// a1 and a2 are alleles (node), can only be HOMOZ_A or HOMOZ_B
bool FounderAlleleGraph4::legal(enum unphased_genotype obs, enum unphased_genotype a1, enum unphased_genotype a2) {
    switch(obs) {
        case HETERO:
            return a1 != a2;
        case HOMOZ_A:
            return (a1 == HOMOZ_A) and (a2 == HOMOZ_A);
        case HOMOZ_B:
            return (a1 == HOMOZ_B) and (a2 == HOMOZ_B);
        case UNTYPED:
            break;
    }
    abort();
}

enum unphased_genotype FounderAlleleGraph4::get_other_allele(enum unphased_genotype obs, enum unphased_genotype a1) {
    switch(obs) {
        case HETERO:
            return a1 == HOMOZ_A ? HOMOZ_B : HOMOZ_A;
        case HOMOZ_A:
            return a1 == HOMOZ_A ? HOMOZ_A : UNTYPED;
        case HOMOZ_B:
            return a1 == HOMOZ_A ? UNTYPED : HOMOZ_A;
        case UNTYPED:
            break;
    }
    abort();
}

// merge component2 into component1
// flip means that you should swap the allele assignments for things that
// were in component2
void FounderAlleleGraph4::combine_components(int component1, int component2, bool flip) {
    enum unphased_genotype tmp;
    
    for(unsigned int i = 0; i < founder_alleles; ++i) {
        if(group_membership[i] == component2) {
            group_membership[i] = component1;
            
            if(flip) {
                tmp = allele_assignment[0][i];
                allele_assignment[0][i] = allele_assignment[1][i];
                allele_assignment[1][i] = tmp;
            }
        }
    }
}

double FounderAlleleGraph4::update_likelihood(unsigned int personid, enum parentage allele) {
    Person* p = ped->get_by_index(personid);
    int old_fa = edge_list[(personid * 2) + allele];
    int tmp = p->get_parentid(allele) * 2;
    int new_fa = ((edge_list[tmp] == old_fa) ? edge_list[tmp + 1] : edge_list[tmp]);
    
    propagate_fa_update(p, old_fa, new_fa);
    
    return 0.0;
}

void FounderAlleleGraph4::propagate_fa_update(Person* p, int old_fa, int new_fa) {
    int tmp = p->get_internalid() * 2;
    
    if((edge_list[tmp] != old_fa) and (edge_list[++tmp] != old_fa)) {
        return;
    }
    
    for(unsigned i = 0; i < p->num_children(); ++i) {
        propagate_fa_update(p->get_child(i), old_fa, new_fa);
    }
    
    edge_list[tmp] = new_fa;
}

