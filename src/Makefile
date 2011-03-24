CXX	:= g++
CXXFLAGS := -g -O0 -Wall -Werror -Wextra -pedantic -std=c++98 -pipe 
#CXXFLAGS := -g -O0 -Wall -Wextra -pedantic -std=c++98 -pipe
# -Weffc++ # too hardcore...
INCLUDES := -I.
LIBS := 

OBJECTS	:= \
    pedigree_parser.o \
    linkage_parser.o \
    map_parser.o \
    simwalk_descent_graph.o \
    founder_allele_graph.o \
    descent_graph.o \
    elimination.o \
    pedigree.o \
    person.o \
    genotype.o \
    genetic_map.o \
    disease_model.o \
    peel_sequence_generator.o \
    rfunction.o \
    simulated_annealing.o \
    markov_chain.o \
    program.o \
    linkage_program.o \
    progress.o \
    peeler.o \
    lod_calculator.o \
    peel_matrix.o \
    main.o

app: $(OBJECTS)
	$(CXX) -o $@ $(LIBS) $(OBJECTS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) app

