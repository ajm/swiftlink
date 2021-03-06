#include "gpu_rfunction.h"


/* the largest matrix possible is 16-dimensions pre-sum because each
 * dimension is always length 4 (2-bits) and the index is a 32-bit int */
int gpu_offsets[] = {
    1 <<  0, 
    1 <<  2, 
    1 <<  4, 
    1 <<  6, 
    1 <<  8, 
    1 << 10, 
    1 << 12,
    1 << 14,
    1 << 16,
    1 << 18,
    1 << 20,
    1 << 22,
    1 << 24,
    1 << 26,
    1 << 28,
    1 << 30
};


// I am going to assume that the length of 'assignment' is to the number of 
// pedigree members and that everything that is not assigned is -1
//
int rfunction_index(struct rfunction* rf, int* assignment, int length) {
    int i, tmp = 0;
    
    for(i = 0; i < (rf->cutset_length - 1); ++i) {
        tmp += (assignment[rf->cutset[i]] * gpu_offsets[i]);
    }
    
    return tmp;
}

// this works out the index for the 'presum_matrix', though with this
// scheme I can always just do:
//      rf->matrix[index - ((0..3) * offset[rf->cutset_length - 1])]
// when I sum over one of the dimensions so long as I ensure that
// the peel_node is always in cutset[cutset_length - 1]
int rfunction_presum_index(struct rfunction* rf, int* assignment, int length) {
    int i, tmp = 0;
    
    for(i = 0; i < rf->cutset_length; ++i) {
        tmp += (assignment[rf->cutset[i]] * gpu_offsets[i]);
    }
    
    return tmp;
}

// from an index, construct the assignment of genotypes to cutset nodes
// 
void rfunction_presum_assignment(struct rfunction* rf, int ind, int* assignment, int length) {
    int index = ind;
    int i;
    
    for(i = 0; i < length; ++i) {
        assignment[i] = -1;
    }
    
    for(i = (rf->cutset_length - 1); i > -1; --i) {
        assignment[rf->cutset[i]] = index / gpu_offsets[i];
        index %= gpu_offsets[i];
    }
}

float rfunction_get(struct rfunction* rf, int* assignment, int length) {
    return rf->matrix[rfunction_index(rf, assignment, length)];
}

float rfunction_trait_prob(struct gpu_state* state, int id, int value, int locus) {
    struct person* p = GET_PERSON(state, id); //pedigree[id];
    
    if(PERSON_ISTYPED(p)) {
        switch(PERSON_GENOTYPE(p, locus)) {
            case GPU_GENOTYPE_AB :
                return ((value == GPU_TRAIT_AB) || (value == GPU_TRAIT_BA)) ? 0.5 : 0.0;
            case GPU_GENOTYPE_AA :
                return (value == GPU_TRAIT_AA) ? 1.0 : 0.0;
            case GPU_GENOTYPE_BB :
                return (value == GPU_TRAIT_BB) ? 1.0 : 0.0;
            default :
                break;
        }
    }
    
    if(! PERSON_ISFOUNDER(p))
        return 0.25;
    
    return MAP_PROB(GET_MAP(state), locus, value);
}

float rfunction_trans_prob(struct gpu_state* state, int locus, int peelnode, 
                           int parent_trait, int child_trait, int parent) {
    
    int trait = get_trait(child_trait, parent);
    int meiosis = 0;
    float tmp = 1.0;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct geneticmap* map = GET_MAP(state);
    
    switch(parent_trait) {
        case GPU_TRAIT_AA:
            return (trait == GPU_TRAIT_A) ? 0.25 : 0.0;
        case GPU_TRAIT_BB:
            return (trait == GPU_TRAIT_B) ? 0.25 : 0.0;
        case GPU_TRAIT_AB:
            meiosis = (trait == GPU_TRAIT_A) ? 0 : 1;
            break;
        case GPU_TRAIT_BA:
            meiosis = (trait == GPU_TRAIT_A) ? 1 : 0;
            break;
    }
    
    if(locus != 0) {
        tmp *= (DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus - 1, parent)) == meiosis ? \
                MAP_INVERSETHETA(map, locus - 1) : \
                MAP_THETA(map, locus - 1));
    }
    
    if(locus != (MAP_LENGTH(map) - 1)) {
        tmp *= (DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, parent)) == meiosis ? \
                MAP_INVERSETHETA(map, locus) : \
                MAP_THETA(map, locus));
    }
    
    return tmp;
}

int get_trait(int value, int parent) {
    switch(parent) {
        case GPU_MATERNAL_ALLELE:
            return (((value == GPU_TRAIT_AA) || (value == GPU_TRAIT_AB)) ? GPU_TRAIT_A : GPU_TRAIT_B);    
        case GPU_PATERNAL_ALLELE:
            return (((value == GPU_TRAIT_AA) || (value == GPU_TRAIT_BA)) ? GPU_TRAIT_A : GPU_TRAIT_B);            
        default:
            break;
    }
    abort();
}

float get_random() {
    return random() / (float) RAND_MAX; 
}

void rfunction_sample(struct rfunction* rf, int* assignment, int assignment_length) {
    float prob_dist[NUM_ALLELES];
    float total = 0.0;
    float r = get_random();
    int peelnode = RFUNCTION_PEELNODE(rf);
    int i;
    
    // extract probabilities
    for(i = 0; i < NUM_ALLELES; ++i) {
        assignment[peelnode] = i;
        prob_dist[i] = RFUNCTION_PRESUM_GET(rf, rfunction_presum_index(rf, assignment, assignment_length));
        total += prob_dist[i];        
    }
    
    // normalise
    for(i = 0; i < NUM_ALLELES; ++i) {
        prob_dist[i] /= total;
    }
    
    /*
    printf("node = %d\n", peelnode);
    printf("random = %f\n", r);
    for(i = 0; i < NUM_ALLELES; ++i) {
        printf("  sample[%d] = %f\n", i, prob_dist[i]);
    }
    */

    // sample
    total = 0.0;
    for(i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        if(r < total) {
            assignment[peelnode] = i;
            return;
        }
    }
    
    abort();
}

void rfunction_evaluate_partner_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    int assignment_length = state->pedigree_length;
    int assignment[assignment_length];
    int peelnode;
    int peelnode_value;
    
    rfunction_presum_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    
    RFUNCTION_PRESUM_SET(rf, ind, \
        rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
        (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
        (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)));

    //printf(" - %d: %e\n", ind, RFUNCTION_PRESUM_GET(rf, ind));
    
    /*
    if(peelnode == 3) {
        printf("\t%f %f %f\n", rfunction_trait_prob(state, peelnode, peelnode_value, locus),
        (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)),
        (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)));
    }
    */
}

void rfunction_evaluate_child_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[assignment_length];
    int peelnode;
    int peelnode_value;
    int mother_value;
    int father_value;
    
    rfunction_presum_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    p = GET_PERSON(state, peelnode);
    mother_value = assignment[PERSON_MOTHER(p)];
    father_value = assignment[PERSON_FATHER(p)];

    RFUNCTION_PRESUM_SET(rf, ind, \
        rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
        rfunction_trans_prob(state, locus, peelnode, mother_value, peelnode_value, GPU_MATERNAL_ALLELE) * \
        rfunction_trans_prob(state, locus, peelnode, father_value, peelnode_value, GPU_PATERNAL_ALLELE) * \
        (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
        (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)));
        
    //printf(" - %d: %e\n", ind, RFUNCTION_PRESUM_GET(rf, ind));
}

void rfunction_evaluate_parent_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[assignment_length];
    int peelnode;
    int peelnode_value;
    int i;
    float tmp;
    
    rfunction_presum_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    
    tmp = rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
          (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
          (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length));
    
    for(i = 0; i < rf->cutset_length; ++i) {
        p = GET_PERSON(state, rf->cutset[i]);
        
        if(! PERSON_ISPARENT(p, peelnode))
            continue;
            
        tmp *= (rfunction_trans_prob(state, locus, PERSON_ID(p), assignment[PERSON_MOTHER(p)], assignment[PERSON_ID(p)], GPU_MATERNAL_ALLELE) * \
                rfunction_trans_prob(state, locus, PERSON_ID(p), assignment[PERSON_FATHER(p)], assignment[PERSON_ID(p)], GPU_PATERNAL_ALLELE));
    }
    
    RFUNCTION_PRESUM_SET(rf, ind, tmp);
}

void rfunction_sum(struct rfunction* rf, int ind) {
    int i;
    int tmp = rf->cutset_length - 1;
    
    RFUNCTION_SET(rf, ind, 0.0);
    
    // unroll?
    for(i = 0; i < 4; ++i) {
        RFUNCTION_ADD(rf, ind, RFUNCTION_PRESUM_GET(rf, ind + (gpu_offsets[tmp] * i)));
    }
    
    //printf(" + %d: %e\n", ind, RFUNCTION_GET(rf, ind));
}

void rfunction_evaluate_element(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    switch(RFUNCTION_TYPE(rf)) {
        case GPU_CHILD_PEEL:
            rfunction_evaluate_child_peel(rf, state, locus, ind);
            break;
            
        case GPU_PARTNER_PEEL:
            rfunction_evaluate_partner_peel(rf, state, locus, ind);
            break;
            
        case GPU_PARENT_PEEL:
            rfunction_evaluate_parent_peel(rf, state, locus, ind);
            break;
            
        default:
            abort();
            break;
    }
}

void rfunction_print(struct rfunction* rf) {
    printf("rfunction: peelnode = %d, type = %d\n", RFUNCTION_PEELNODE(rf), RFUNCTION_TYPE(rf));
}

// for now, to test the code on the CPU
void rfunction_evaluate(struct rfunction* rf, struct gpu_state* state, int locus) {
    int i;
    
    //rfunction_print(rf);
    
    for(i = 0; i < rf->presum_length; ++i) {
        rfunction_evaluate_element(rf, state, locus, i);
    }
    
    for(i = 0; i < rf->matrix_length; ++i) {
        rfunction_sum(rf, i);
    }
}

// i don't know how this is actually going to work
// i specify the grid/block/thread layout, and then each thread
// should correspond to a single cell in a matrix, ie: in the
// presum matrix, however this is a problem, because there are 
// a different number of cells in each r-function, i don't know 
// if I just state the max needed for a given problem and there 
// is a way to kill the ones that fall out of the range or 
// whether the thing to do would be to calculate the index modulo
// the number of threads?
//void rfunction_kernel_start(struct gpu_state* state) {
    // figure out what block/thread we are
    // + what locus that relates to
    // a) loop through all r-functions for that locus
    //      i)  calculate presum_matrix
    //      ii) calculate summed matrix
    // b) synchronise all threads
    // c) do the same thing for this locus + 1    
//}

int sample_hetero_mi(int allele, int trait) {
    if(allele == GPU_TRAIT_A) {
        return (trait == GPU_TRAIT_AB) ? 0 : 1;
    }
    else {
        return (trait == GPU_TRAIT_AB) ? 1 : 0;
    }
}

// find prob of setting mi to 0
// find prob of setting mi to 1
// normalise + sample
int sample_homo_mi(struct gpu_state* state, int personid, int locus, int parent) {
    float prob_dist[2];
    float total;
    int i;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct geneticmap* map = GET_MAP(state);
    
    prob_dist[0] = 1.0;
    prob_dist[1] = 1.0;
    
    if(locus != 0) {
        i = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus - 1, parent));
        prob_dist[0] *= (i == 0 ? MAP_INVERSETHETA(map, locus - 1) : MAP_THETA(map, locus - 1));
        prob_dist[1] *= (i == 1 ? MAP_INVERSETHETA(map, locus - 1) : MAP_THETA(map, locus - 1));
    }
    
    // map length is actually the number of recombination fractions (#markers - 1)
    if(locus != (MAP_LENGTH(map) - 1)) {
        i = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus + 1, parent));
        prob_dist[0] *= (i == 0 ? MAP_INVERSETHETA(map, locus) : MAP_THETA(map, locus));
        prob_dist[1] *= (i == 1 ? MAP_INVERSETHETA(map, locus) : MAP_THETA(map, locus));
    }
    
    total = prob_dist[0] + prob_dist[1];
    
    prob_dist[0] /= total;
    prob_dist[1] /= total;
    
    return (get_random() < prob_dist[0]) ? 0 : 1;
}

// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right    
int sample_mi(struct gpu_state* state, int allele, int trait, int personid, int locus, int parent) {
    switch(trait) {
        case GPU_TRAIT_AB:
        case GPU_TRAIT_BA:
            return sample_hetero_mi(allele, trait);
            
        case GPU_TRAIT_AA:
        case GPU_TRAIT_BB:
            return sample_homo_mi(state, personid, locus, parent);
            
        default:
            abort();
    }
}

// sample meiosis indicators
// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right
void sample_meiosis_indicators(struct gpu_state* state, int* assignment, int locus) {
    struct person* p;
    struct descentgraph* dg;
    int i;
    int mother_trait;
    int father_trait;
    int person_trait;
    int mat_allele;
    int pat_allele;
    int mat_mi;
    int pat_mi;
    
    dg = GET_DESCENTGRAPH(state);
    
    for(i = 0; i < state->pedigree_length; ++i) {
        p = GET_PERSON(state, i);
        
        if(PERSON_ISFOUNDER(p))
            continue;
        
        person_trait = assignment[p->id];
        mother_trait = assignment[p->mother];
        father_trait = assignment[p->father];
        
        mat_allele = (person_trait == GPU_TRAIT_AA) || (person_trait == GPU_TRAIT_AB) ? GPU_TRAIT_A : GPU_TRAIT_B;
        pat_allele = (person_trait == GPU_TRAIT_AA) || (person_trait == GPU_TRAIT_BA) ? GPU_TRAIT_A : GPU_TRAIT_B;
        
        mat_mi = sample_mi(state, mat_allele, mother_trait, i, locus, GPU_MATERNAL_ALLELE);
        pat_mi = sample_mi(state, pat_allele, father_trait, i, locus, GPU_PATERNAL_ALLELE);
        
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, i, locus, GPU_MATERNAL_ALLELE), mat_mi);
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, i, locus, GPU_PATERNAL_ALLELE), pat_mi);
    }
}

void sampler_step(struct gpu_state* state, int locus) {
    int i;
    int assignment[state->pedigree_length];
    
    // forward peel
    for(i = 0; i < state->functions_per_locus; ++i) {
        rfunction_evaluate(GET_RFUNCTION(state, i, locus), state, locus);
    }
    
    // reverse peel, sampling ordered genotypes
    for(i = state->functions_per_locus - 1; i >= 0; --i) {
        rfunction_sample(GET_RFUNCTION(state, i, locus), assignment, state->pedigree_length);
    }
    
    sample_meiosis_indicators(state, assignment, locus);
}

void print_ints(int* data, int length) {
    int i;
    
    for(i = 0; i < length; ++i) {
        printf("%d ", data[i]);
    }
    
    printf("\n");
}

void print_rfunction(struct rfunction* r) {
    printf( "\tid: %d\n"
            "\ttype: %d\n"
            "\tprev1: %d\n"
            "\tprev2: %d\n",
            r->id,
            r->peel_type,
            r->prev1 != NULL ? r->prev1->id : -1,
            r->prev2 != NULL ? r->prev2->id : -1
    );
    
    printf("\tcutset: ");
    print_ints(r->cutset, r->cutset_length);
    printf("\n");
}

void print_person(struct person* p) {
    printf( "\tid: %d\n"
            "\tmother: %d\n"
            "\tfather: %d\n"
            "\tfounder: %d\n"
            "\ttyped: %d\n"
            "\tprobs: %.3f %.3f %.3f %.3f\n",
            p->id,
            p->mother,
            p->father,
            p->isfounder,
            p->istyped,
            p->prob[0],
            p->prob[1],
            p->prob[2],
            p->prob[3]
    );
    
    printf("\tgenotypes: ");
    print_ints(p->genotypes, p->genotypes_length);
    printf("\n");
}

void print_map(struct geneticmap* map) {
    int i, j;
    
    for(i = 0; i < (map->map_length - 1); ++i) {
        printf("\t%d\t%.3f %.3f\n", i, map->thetas[i], map->inversethetas[i]);
    }
    printf("\n");
    
    for(i = 0; i < map->map_length; ++i) {
        printf("\t%d:\t", i);
        for(j = 0; j < 4; ++j) {
            printf("%.3f ", MAP_PROB(map, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void print_descentgraph(struct descentgraph* dg, int ped_length, int map_length) {
    int i, j;
    
    for(i = 0; i < ped_length; ++i) {
        printf("\t%d:\t", i);
        for(j = 0; j < map_length; ++j) {
            printf( "%d%d ",
                    DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, j, GPU_MATERNAL_ALLELE)),
                    DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, j, GPU_PATERNAL_ALLELE))
            );
        }
        printf("\n");
    }
    printf("\n");
}

void print_everything(struct gpu_state* state) {
    int i;
    
    printf("RFUNCTIONS:\n");
    for(i = 0; i < state->functions_length; ++i) {
        print_rfunction(&state->functions[i]);
    }
    printf("\n");
    
    printf("PEDIGREE:\n");
    for(i = 0; i < state->pedigree_length; ++i) {
        print_person(&state->pedigree[i]);
    }
    printf("\n");
    
    printf("MAP:\n");
    print_map(state->map);
    printf("\n");
    
    printf("DESCENTGRAPH:\n");
    print_descentgraph(state->dg, state->pedigree_length, state->map->map_length);
    printf("\n");
}

