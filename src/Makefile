CXX	:= g++
CC := gcc
#CXXFLAGS := -g -O3 -Wall -Wextra -pedantic -std=c++98 -pipe -Wno-variadic-macros -Weffc++ -Werror
#CXXFLAGS := -g -O2 -Wall -Wextra -pedantic -std=c++98 -pipe -Wno-variadic-macros -Weffc++
CXXFLAGS := -g -O2 -Wall -Wextra -Weffc++ -pedantic -std=c++98 -pipe
CFLAGS := -g -O2 -Wall -pipe
INCLUDES := -I.
LIBS := 

OBJECTS	:= \
	trait.o \
	genotype.o \
	types.o \
	logarithms.o \
	disease_model.o \
	genetic_map.o \
	person.o \
	pedigree.o \
	descent_graph.o \
	pedigree_parser.o \
	linkage_parser.o \
	map_parser.o \
	linkage_writer.o \
	haplotype_writer.o \
	peel_sequence_generator.o \
	founder_allele_graph.o \
	founder_allele_graph2.o \
	elimination.o \
	peel_matrix.o \
	rfunction.o \
	trait_rfunction.o \
	sampler_rfunction.o \
	peeler.o \
	lod_calculator.o \
	markov_chain.o \
	locus_sampler2.o \
	meiosis_sampler.o \
	linkage_program.o \
	program.o \
	progress.o \
	main.o \
	gpu_wrapper.o \
	gpu_rfunction.o \
	test_program.o

app: $(OBJECTS)
	$(CXX) -o $@ $(LIBS) $(OBJECTS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) app

