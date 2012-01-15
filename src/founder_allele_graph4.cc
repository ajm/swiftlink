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
    Person* p;
    int tmp, pid;
    
    for(unsigned int i = 0; i < ped->num_members(); ++i) {
	    pid = (*sequence)[i];
        p = ped->get_by_index(pid);
	    tmp = pid * 2;
	    
	    ss << pid << "\t(" << edge_list[tmp] << ", " << edge_list[tmp + 1] << ")\n";
	}
    
    return ss.str();
}

double FounderAlleleGraph4::likelihood() {
    Person* p;
    int tmp;
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
    group_active.assign(founder_alleles, false);
    
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
        if(mat_fa == pat_fa) {
            if(g == HETERO) {
                return 0.0;
            }
            
            group1 = group_membership[mat_fa];
            
            if(group1 != DEFAULT_COMPONENT) {
                // already belongs to a component
                fixed1 = group_fixed[group1];
                if(fixed1 != -1) {
                    if(g != allele_assignment[fixed1][mat_fa]) {
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
                        return 0.0;
                    }
                }
            }
            else {
                group_membership[mat_fa] = group_index;
                group_fixed[group_index] = 0;
                group_active[group_index] = true;
                group_size[group_index] = 1;
                allele_assignment[0][mat_fa] = g;
                prob[0][group_index] = (g == HOMOZ_A ? major_freq : minor_freq);
                
                ++group_index;
                ++num_groups;
            }
        }
        // not autozygous
        else {
            group1 = group_membership[mat_fa];
            group2 = group_membership[pat_fa];
            
            if(group1 != DEFAULT_COMPONENT) {
                
                fixed1 = group_fixed[group1];
                
                if(group2 != DEFAULT_COMPONENT) {
                    // both in the same group
                    if(group1 == group2) {
                        if(fixed1 != -1) {
                            // fixed, check if legit
                            if(not legal(g, allele_assignment[fixed1][mat_fa], allele_assignment[fixed1][pat_fa])) {
                                return 0.0;
                            }
                        }
                        else {
                            // not fixed, check if still the case
                            legal0 = legal(g, allele_assignment[0][mat_fa], allele_assignment[0][pat_fa]);
                            legal1 = legal(g, allele_assignment[1][mat_fa], allele_assignment[1][pat_fa]);
                            
                            if(legal0) {
                                if(not legal1) {
                                    group_fixed[group1] = 0;
                                }
                            }
                            else {
                                if(legal1) {
                                    group_fixed[group1] = 1;
                                }
                                else {
                                    return 0.0;
                                }
                            }
                        }
                    }
                    else {
                        // in different groups
                        fixed2 = group_fixed[group2];
                        
                        if(fixed1 != -1) {
                            if(fixed2 != -1) {
                                // fixed, check if legit
                                if(not legal(g, allele_assignment[fixed1][mat_fa], allele_assignment[fixed2][pat_fa])) {
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
                            else {
                                // all true, assignments don't need swapping
                                fixed1 = fixed2 = -1;
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
                            return 0.0;
                        }
                        else {
                            allele_assignment[fixed1][pat_fa] = tmp;
                            prob[fixed1][group1] *= (tmp == HOMOZ_A ? major_freq : minor_freq);
                        }
                    }
                    else {
                        enum unphased_genotype tmp0 = get_other_allele(g, allele_assignment[0][mat_fa]);
                        enum unphased_genotype tmp1 = get_other_allele(g, allele_assignment[1][mat_fa]);
                        
                        if(tmp0 != UNTYPED) {
                            if(tmp1 != UNTYPED) {
                                allele_assignment[0][pat_fa] = tmp0;
                                allele_assignment[1][pat_fa] = tmp1;
                                prob[0][group1] *= (tmp0 == HOMOZ_A ? major_freq : minor_freq);
                                prob[1][group1] *= (tmp1 == HOMOZ_A ? major_freq : minor_freq);
                            }
                            else {
                                allele_assignment[0][pat_fa] = tmp0;
                                prob[0][group1] *= (tmp0 == HOMOZ_A ? major_freq : minor_freq);
                                group_fixed[group1] = 0;
                            }
                        }
                        else {
                            if(tmp1 != UNTYPED) {
                                allele_assignment[1][pat_fa] = tmp1;
                                prob[1][group1] *= (tmp1 == HOMOZ_A ? major_freq : minor_freq);
                                group_fixed[group1] = 1;
                            }
                            else {
                                return 0.0;
                            }
                        }
                    }
                    
                    group_membership[pat_fa] = group1;
                    group_size[group1] += 1;
                }
            }
            else if(group2 != DEFAULT_COMPONENT) {
                fixed2 = group_fixed[group2];
                
                if(fixed2 != -1) {
                    enum unphased_genotype tmp = get_other_allele(g, allele_assignment[fixed2][pat_fa]);
                    
                    if(tmp == UNTYPED) {
                        return 0.0;
                    }
                    else {
                        allele_assignment[fixed2][mat_fa] = tmp;
                        prob[fixed2][group2] *= (tmp == HOMOZ_A ? major_freq : minor_freq);
                    }
                }
                else {
                    enum unphased_genotype tmp0 = get_other_allele(g, allele_assignment[0][pat_fa]);
                    enum unphased_genotype tmp1 = get_other_allele(g, allele_assignment[1][pat_fa]);
                    
                    if(tmp0 != UNTYPED) {
                        if(tmp1 != UNTYPED) {
                            allele_assignment[0][mat_fa] = tmp0;
                            allele_assignment[1][mat_fa] = tmp1;
                            prob[0][group2] *= (tmp0 == HOMOZ_A ? major_freq : minor_freq);
                            prob[1][group2] *= (tmp1 == HOMOZ_A ? major_freq : minor_freq);
                        }
                        else {
                            allele_assignment[0][mat_fa] = tmp0;
                            prob[0][group2] *= (tmp0 == HOMOZ_A ? major_freq : minor_freq);
                            group_fixed[group2] = 0;
                        }
                    }
                    else {
                        if(tmp1 != UNTYPED) {
                            allele_assignment[1][mat_fa] = tmp1;
                            prob[1][group2] *= (tmp1 == HOMOZ_A ? major_freq : minor_freq);
                            group_fixed[group2] = 1;
                        }
                        else {
                            return 0.0;
                        }
                    }
                }
                
                group_membership[mat_fa] = group2;
                group_size[group2] += 1;
            }
            else {
                if(g == HETERO) {
                    allele_assignment[0][mat_fa] = allele_assignment[1][pat_fa] = HOMOZ_A;
                    allele_assignment[1][mat_fa] = allele_assignment[0][pat_fa] = HOMOZ_B;
                    prob[0][group_index] = \
                    prob[1][group_index] = major_freq * minor_freq;
                    group_fixed[group_index] = -1;
                }
                else {
                    allele_assignment[0][mat_fa] = allele_assignment[0][pat_fa] = g;
                    prob[0][group_index] = get_freq(g) * get_freq(g);
                    group_fixed[group_index] = 0;
                }
                
                // neither in a group
                group_membership[mat_fa] = \
                group_membership[pat_fa] = group_index;
                group_size[group_index] = 2;
                group_active[group_index] = true;
                ++group_index;
                ++num_groups;
            }
        }
    }
    
    // calculate likelihood
    double ret_prob = 1.0;
    
    for(int i = 0; i < group_index; ++i) {
        if(group_active[i] and (group_size[i] > 1)) {
            fixed1 = group_fixed[i];
            if(fixed1 != -1) {
                ret_prob *= prob[fixed1][i];
            }
            else {
                ret_prob *= (prob[0][i] + prob[1][i]);
            }
        }
    }
    
    return ret_prob;
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
            return a1 == HOMOZ_B ? HOMOZ_B : UNTYPED;
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
    
    if(flip) {
        prob[0][component1] *= prob[1][component2];
        prob[1][component1] *= prob[0][component2];
    }
    else {
        prob[0][component1] *= prob[0][component2];
        prob[1][component1] *= prob[1][component2];
    }
    
    group_size[component1] += group_size[component2];
    
    group_active[component2] = false;
}

void FounderAlleleGraph4::reset(DescentGraph& dg) {
    Person* p;
    int pid;
    int parent_allele;
    int tmp;
    
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
}

void FounderAlleleGraph4::flip(unsigned int personid, enum parentage allele) {
    Person* p = ped->get_by_index(personid);
    int old_fa = edge_list[(personid * 2) + allele];
    int tmp = p->get_parentid(allele) * 2;
    int new_fa = ((edge_list[tmp] == old_fa) ? edge_list[tmp + 1] : edge_list[tmp]);
    
    if(old_fa == new_fa)
        return;
    
    propagate_fa_update(p, allele, old_fa, new_fa);
}

void FounderAlleleGraph4::propagate_fa_update(Person* p, enum parentage allele, int old_fa, int new_fa) {
    int tmp = p->get_internalid() * 2;
    
    if(edge_list[tmp + allele] != old_fa)
        return;
    
    enum parentage newallele = p->isfemale() ? MATERNAL : PATERNAL;
    
    for(unsigned i = 0; i < p->num_children(); ++i) {
        propagate_fa_update(p->get_child(i), newallele, old_fa, new_fa);
    }
    
    edge_list[tmp + allele] = new_fa;
}

