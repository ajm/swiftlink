CXX	:= g++
#CXXFLAGS := -g -O0 -Wall -Werror -Wextra -pedantic -std=c++98 -pipe 
CXXFLAGS := -g -O3 -Wall -Wextra -pedantic -std=c++98 -pipe
# -Weffc++ # too hardcore...
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
	peeler.o \
	peel_matrix.o \
	lod_calculator.o \
	simwalk_sampler.o \
	elimination.o \
	progress.o \
	simulated_annealing.o \
    markov_chain.o \
    linkage_program.o \
	program.o \
	main.o

app: $(OBJECTS)
	$(CXX) -o $@ $(LIBS) $(OBJECTS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) app

