CXX	:= g++
CC := gcc
GPUCC := nvcc

# 32bit
#CXXFLAGS := -g -O2 -Wall -Wextra -Weffc++ -pedantic -std=c++98 -pipe -fopenmp -m32
#CFLAGS := -g -O2 -Wall -pipe -fopenmp -m32
#LDFLAGS := -m32
#INCLUDES := -I. -I/usr/local/cuda/include
#LIBS := -L/usr/local/cuda/lib -L/usr/lib32 -lcudart -lgomp -lgsl -lgslcblas -lrt
#GPUFLAGS := -arch=sm_20 -O2 -m32

CXXFLAGS := -g -O2 -Wall -Wextra -Weffc++ -pedantic -std=c++98 -pipe -fopenmp -DUSE_CUDA #-DCUDA_SHAREDMEM_CACHE
CFLAGS := -g -O2 -Wall -pipe -fopenmp
LDFLAGS :=
INCLUDES := -I. -I/usr/local/cuda/include
LIBS := -lgomp -lgsl -lgslcblas -lrt -L/usr/local/cuda/lib64 -lcudart
GPUFLAGS := -arch=sm_20 -O2

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
	gpu_lodscores.o \
	cuda_linkage.o

swift: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cu
	$(GPUCC) $(GPUFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) swift

# Dependencies for File: SwiftLink:
SwiftLink: "Allegro "Linkage (anecdotally, ), /github.com/ajm/swiftlink.git
SwiftLink: /stackoverflow.com/questions/10320562/a-disastrous-slowdown-of-cudamalloc-in-nvidia-drivers-from-version-285
SwiftLink: 32-bit 4.0 5.0 </b> Allegro. CUDA Format" I If LINKAGE MCMC Mega2
SwiftLink: Parallel RC1 We
SwiftLink: [Alohomora](http:/bioinformatics.oxfordjournals.org/content/21/9/2123.full.pdf)
SwiftLink: analysis and are bug). but by compatible confirms data different.
SwiftLink: does example file file, files files, files. fine for format, format.
SwiftLink: from generate generated generates have however in input is it
SwiftLink: linkage locus map mode mostly must not on pedigree seem seems select
SwiftLink: selecting simple something standard suffer tested that the this to
SwiftLink: using very work written you your
# Dependencies for File: <h1>SwiftLink:
<h1>SwiftLink: (anecdotally, /github.com/ajm/swiftlink.git
<h1>SwiftLink: /stackoverflow.com/questions/10320562/a-disastrous-slowdown-of-cudamalloc-in-nvidia-drivers-from-version-285">slow
<h1>SwiftLink: 32-bit 4.0 5.0 </p> <a Allegro.</p> CUDA LINKAGE MCMC Parallel
<h1>SwiftLink: RC1 We analysis</h1> and bug bug).</p> by cudaMalloc</a>, data
<h1>SwiftLink: does down file file, files fine for format. from generated have
<h1>SwiftLink: however
<h1>SwiftLink: href="http:/bioinformatics.oxfordjournals.org/content/21/9/2123.full.pdf">Alohomora</a>
<h1>SwiftLink: identical in input it linkage map mode mostly not on pedigree
<h1>SwiftLink: seem seems standard suffer tested this to using work
