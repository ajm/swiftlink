#!/usr/bin/env python
import sys
import random
import copy

#debug = True
debug = False

class PeelOperation :
    def __init__(self, s, pivot, cutset) :
        self.s = s
        self.pivot = pivot
        self.cutset = cutset

    def __str__(self) :
        return "%s %s cutset=(%s)" % (self.s, str(self.pivot), ', '.join(map(str, self.cutset)))

class R_function :
    def __init__(self, pivot, cutset, node) :
        self.pivot = pivot
        self.cutset = cutset
        self.matrix = {} # size = 3**len(cutset)
        self.node = node # node that has this R-function as a term, 
                         # ie: next node dealt with in cutset
    
    def __str__(self) :
        pass

class PedInputError(Exception) :
    pass

class Person :
    # allele
    UNDEFINED_ALLELE = 0
    # sex
    UNDEFINED = 0
    MALE = 1
    FEMALE = 2
    # affection status
    UNDEFINED_AFFECTION = 0
    UNAFFECTED = 1
    AFFECTED = 2

    FOUNDER = 0

    WHITE = 0
    GREY = 1
    BLACK = 2

    HOMO_A = 0
    HETERO = 1
    HOMO_B = 2
    
    
    def __init__(self, id, father, mother, sex, affected, ped) :
        self.id = id        
        self.father = father
        self.mother = mother
        self.sex = sex
        self.affected = affected

        self.ped = ped
        
        self.offspring = []
        self.partners = []

        self.genotypes = []

        self.peeled = False
        
        if (self.mother == 0) ^ (self.father == 0) :
            raise PedInputError("one parent was undefined")

    def add_partner(self, id) :
        if id not in self.partners :
            self.partners.append(id)

    def add_offspring(self, id) :
        if id not in self.offspring :
            self.offspring.append(id)

    def has_offspring(self) :
        return len(self.offspring) != 0

    def is_affected(self) :
        return self.affected == Person.AFFECTED

    def is_leaf(self) :
        return not self.has_offspring()

    def is_founder(self) :
        return self.father == 0 and self.mother == 0

    def get_other_parent(self, id) :
        if id == self.mother :
            return self.father
        elif id == self.father :
            return self.mother
        else :
            raise Exception('bad parent request')

    def _peeled(self, people) :
        return map(lambda x : self.ped[x].peeled, people).count(False)

    def offspring_unpeeled(self) :
        return self._peeled(self.offspring)

    def parents_unpeeled(self) :
        if self.is_founder() :
            return 0

        return self._peeled([self.father, self.mother])

    def partners_unpeeled(self) :
        return self._peeled(self.partners)

    # both parents have been peeled
    def ripe_above(self) :
        return self.parents_unpeeled() == 0

    # all children have been peeled, if they exist
    # partner has been peeled, if they exist
    def ripe_below(self) :
        #return (self.offspring_unpeeled() == 0) and (self.partners_unpeeled() == 0)
        return self.offspring_unpeeled() == 0

    def ripe_across(self) :
        return (self.offspring_unpeeled() == 0) and (self.parents_unpeeled() == 0) and (self.partners_unpeeled() == 1)

    def ripe_all(self) :
        return (self.offspring_unpeeled() == 0) and (self.parents_unpeeled() == 0) and (self.partners_unpeeled() == 0)

    def at_least_one_parent_peripheral(self) :
        return True in map(lambda x : self.ped[x].ripe_above(), [self.mother, self.father])

    def generic_peel(self) :
        # last person, just sum
        if self.ripe_all() :
            return PeelOperation("Last", self.id, self.get_cutset())
                #[])

        # peel partners
        if self.ripe_across() :
            return PeelOperation("PeelPartner", self.id, self.get_cutset()) #[self.partners[0]])

        # peel a child on to parents and unpeeled partners
        # XXX this is not correct if this person is not related to their partner
        # somehow, however this will work because there will always be a peripheral
        # family with a smaller cutset if there is not an inbreeding loop
        if self.ripe_below() :
            return PeelOperation("PeelUp", self.id, self.get_cutset())
                #[self.father, self.mother] + filter(lambda x : not self.ped[x].peeled, self.partners))

        # peel parents on to all unpeeled children
        if self.ripe_above() and self.ped[self.partners[0]].ripe_above() :
            return PeelOperation("PeelDown", self.id, self.get_cutset())
                #[self.id, self.partners[0]] + filter(lambda x : not ped[x].peeled, self.offspring))
        
        # return empty list, 
        # (XXX or maybe there are multiple possibilities? esp with a loop)
        return None

    def neighbours(self) :
        return filter(lambda x : x != Person.FOUNDER, [self.mother, self.father] + self.partners + self.offspring)

    def get_cutset(self) :
        # perform bfs, stopping where there are unpeeled people
        # store these people in a list and return
        stack = [self.id]
        visited = {}
        cutset = []

        for i in self.ped :
            visited[i] = Person.WHITE
        visited[self.id] = Person.GREY
        
        while len(stack) != 0 :
            tmp = stack.pop()
            for i in self.ped[tmp].neighbours() :
                if not self.ped[i].peeled :
                    if (i not in cutset) and (i != self.id) :
                        cutset.append(i)
                    continue

                if visited[i] == Person.WHITE :
                    stack.append(i)
                    visited[i] = Person.GREY
            
            visited[tmp] = Person.BLACK

        return cutset


def read_pedfile(fname) :
    f = open(fname)
    p = {}

    for line in f :
        line = line.strip()
        if len(line) == 0 :
            continue
        data = line.split()
        info = map(int, data[:6])
        p[info[1]] = Person(info[1], info[2], info[3], info[4], info[5], p)
        
        # read a single genotype
        if (info[6] == '0') or (info[7] == '0') :
            p[info[1]].genotypes = [0.25, 0.5, 0.25]
        elif info[6] != info[7] :
            p[info[1]].genotypes = [0.0, 1.0, 0.0]
        else :
            if p[info[6] == '1' :
                p[info[1]].genotypes = [1.0, 0.0, 0.0]
            else :
                p[info[1]].genotypes = [0.0, 0.0, 1.0]
        
    f.close()
    
    return p

def fill_in_details(ped) :
    for k in ped :
        p = ped[k]
        if not p.is_founder() :
            ped[p.mother].add_offspring(k)
            ped[p.father].add_offspring(k)
            ped[p.mother].add_partner(p.father)
            ped[p.father].add_partner(p.mother)
    
    return ped

def create_peelorder(ped) :
    peelorder = []

    while True :
        tmp = []
        unpeeled_count = 0
        for k in ped :
            if not ped[k].peeled : 
                unpeeled_count += 1
                tmp2 = ped[k].generic_peel()
                if tmp2 != None :
                    tmp.append(tmp2)
        
        
        if unpeeled_count == 0 :
            break
        
        
        if len(tmp) != 0 :
            cutset_size = min(map(lambda x : len(x.cutset), tmp))            
            tmp2 = filter(lambda x : len(x.cutset) == cutset_size, tmp)
            
            random.shuffle(tmp2) # XXX <--- just to see if robust
            
            if debug :
                for i in tmp :
                    print "\t%s" % str(i)
                print "\n* %s\n" % str(tmp2[0])
            
            peelorder.append(tmp2[0])
            ped[tmp2[0].pivot].peeled = True
            continue

    return peelorder

def print_peelorder(peelorder) :
    for i in peelorder :
        print str(i)

if __name__ == '__main__' :
    ped = fill_in_details(read_pedfile(sys.argv[1]))
    peelorder = create_peelorder(ped)
    print_peelorder(peelorder)


