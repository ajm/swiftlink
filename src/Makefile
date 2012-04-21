CXX	:= g++
CC := gcc
GPUCC := nvcc
CXXFLAGS := -g -O3 -Wall -Wextra -Weffc++ -pedantic -std=c++98 -pipe -fopenmp -m32
CFLAGS := -g -O3 -Wall -pipe -fopenmp -m32
LDFLAGS := -m32
INCLUDES := -I. -I/usr/local/cuda/include
#LIBS := -lgomp -lgsl -lgslcblas -L/usr/local/cuda/lib64 -lcudart
LIBS := -L/usr/local/cuda/lib -L/usr/lib32 -lcudart -lgomp -lgsl -lgslcblas -lrt
GPUFLAGS := -arch=sm_20 -O3 -m32

OBJECTS	:= \
	trait.o \
	genotype.o \
	types.o \
	random.o \
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
	peel_sequence_generator.o \
	founder_allele_graph4.o \
	founder_allele_graph3.o \
	elimination.o \
	peel_matrix.o \
	rfunction.o \
	sampler_rfunction.o \
	trait_rfunction.o \
	peeler.o \
	markov_chain.o \
	sequential_imputation.o \
	locus_sampler2.o \
	meiosis_sampler.o \
	linkage_program.o \
	program.o \
	progress.o \
	main.o \
	gpu_markov_chain.o \
	cuda_linkage.o \
	tinymt/tinymt32_host.o

apptest: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cu
	$(GPUCC) $(GPUFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) #app

