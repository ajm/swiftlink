CXX	:= g++
#CXXFLAGS := -g -O3 -Wall -Wextra -pedantic -std=c++98 -pipe -Wno-variadic-macros -Weffc++ -Werror
CXXFLAGS := -g -O3 -Wall -Wextra -pedantic -std=c++98 -pipe -Wno-variadic-macros -Weffc++
INCLUDES := -I.
LIBS := 

OBJECTS	:= \
    pedigree_parser.o \
    linkage_parser.o \
    map_parser.o \
	linkage_writer.o \
	haplotype_writer.o \
	genotype.o \
	person.o \
	pedigree.o \
	genetic_map.o \
	disease_model.o \
	descent_graph.o \
    founder_allele_graph.o \
    peel_sequence_generator.o \
    rfunction.o \
	trait_rfunction.o \
	sampler_rfunction.o \
	peeler.o \
	peel_matrix.o \
	lod_calculator.o \
	simwalk_sampler.o \
	elimination.o \
	progress.o \
	simulated_annealing.o \
    markov_chain.o \
	locus_sampler.o \
    linkage_program.o \
	program.o \
	main.o

app: $(OBJECTS)
	$(CXX) -o $@ $(LIBS) $(OBJECTS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) app

