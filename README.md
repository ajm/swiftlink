# SwiftLink: Parallel MCMC linkage analysis

SwiftLink performs multipoint parametric linkage analysis on large consanguineous pedigrees and is primarily targeted at pedigrees that cannot be analysed by a Lander-Green algorithm based program, i.e. many markers, but larger pedigrees. The current version of SwiftLink only supports SNP markers.

The SwiftLink source code is licensed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

If you use SwiftLink in your work please cite [Medlar et al, 2013](http://bioinformatics.oxfordjournals.org/content/29/4/413.long)

    @article{medlar2013swiftlink,
      title={SwiftLink: parallel MCMC linkage analysis using multicore CPU and GPU},
      author={Medlar, Alan and G{\l}owacka, Dorota and Stanescu, Horia and Bryson, Kevin and Kleta, Robert},
      journal={Bioinformatics},
      volume={29},
      number={4},
      pages={413--419},
      year={2013},
      publisher={Oxford Univ Press}
    }

## Installation

SwiftLink's only mandatory dependency is GNU scientific library for the Mersenne Twister pseudo random number generator. Optionally, SwiftLink can be compiled with CUDA support.

Download source code:

    git clone git://github.com/ajm/swiftlink.git

Build without CUDA support:

    cd swiftlink/src
    make

Build with CUDA support:

    cd swiftlink/src
    make -f Makefile.cuda

Build under Mac OS (using [homebrew](http://brew.sh/) for dependencies):

    brew install gsl libiomp clang-omp
    cd swiftlink/src
    make -f Makefile.macos

## Input Files

SwiftLink expects three input files: pedigree file, map file and locus data file, in LINKAGE format. We have mostly tested it on input files generated by [Alohomora](http://bioinformatics.oxfordjournals.org/content/21/9/2123.full.pdf) for Allegro.

> <b>Update:</b> If you are using Mega2 to generate your input files, you must select "Allegro Format" to generate compatible pedigree and locus data files. 
> Mega2 seems to confirm that the examples I have written (found in the example directory) when used as input is in LINKAGE format, but selecting "Linkage Format" as the output format generates something incompatible.

## CUDA versions

The current version of SwiftLink has been tested on Linux with CUDA version 7.5 (tested Sept. 2016).

Older versions of SwiftLink are known to not work properly with CUDA versions 4.1 and 4.2 due to a known [slow down bug in cudaMalloc](http://stackoverflow.com/questions/10320562/a-disastrous-slowdown-of-cudamalloc-in-nvidia-drivers-from-version-285) (anecdotally, CUDA 4.0 in 32-bit mode did not seem to suffer from this bug).

## Examples

All the input files used in the following commands can be found in the [examples](https://github.com/ajm/swiftlink/tree/master/examples) directory. We use the pedigree from [Bockenhauer et al, 2009](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398803/), but the data and map are simulated.

### Expected LOD (ELOD) score

SwiftLink can calculate an expected LOD (ELOD) score for your pedigree assuming a recessive trait with complete penetrance (the default):

    swift -p east.ped --elod

For an X-linked recessive trait with complete penetrance:

    swift -p xlinked.ped --elod -X

For a dominant trait with complete penetrance (for illustrative purposes only, the file dominant.ped is not provided):

    swift -p dominant.ped --elod --penetrance=0.0,1.0,1.0

The ELOD score can be used both as a power analysis and as an additional quality control post-analysis. If the maximum LOD score differs considerably from the ELOD, then this could point to an unidentified problem with the input data or other model misspecification.

### Linkage analysis

The simplest way to run a linkage analysis with SwiftLink, i.e. with default parameters, is the following:

    swift -p east.ped -m east.map -d east.dat -o results.txt

This will perform either an autosomal or X-linked analysis dependent on whether it is specified in the first line of the DAT file (SwiftLink can be forced to perform an X-linked analysis with the -X flag, see options). By default SwiftLink only uses a single CPU core and only performs a single replicate.

### Using multiple CPUs

SwiftLink can be efficiently run across multiple CPU cores. Here we perform the same analysis using four CPUs:

    swift -p east.ped -m east.map -d east.dat -o results.txt -c 4 

### Performing multiple runs

SwiftLink has a builtin function to run multiple Markov chains and output LOD scores averaged over all replicates. For a majority of projects we have been involved in ~10 replicates is sufficient:

    swift -p east.ped -m east.map -d east.dat -o results.txt -c 4 -R 10

### Affected-only analysis

SwiftLink can easily perform an affected-only analysis, forcing all negative affection statuses to unknown status:

    swift -p east.ped -m east.map -d east.dat -o results.txt -c 4 -a

### Using the GPU

If you have a CUDA-compatible GPU and have the CUDA drivers installed (see [CUDA installation guide](http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-linux/)), SwiftLink can offload LOD score calculations to the GPU and speed up the overall runtime. The GPU code only supports autosomal linkage analysis:

    swift -p east.ped -m east.map -d east.dat -o results.txt -c 4 -g

### MCMC options and convergence diagnostics

The [examples](https://github.com/ajm/swiftlink/tree/master/examples) directory contains mostly toy examples, but depending on the complexity of your project you may need to spend some time ensuring that the Markov chain has converged to the stationary distribution to ensure your inferences are trustworthy. 

This command runs SwiftLink for 1,000,000 iterations of burnin, followed by 1,000,000 iterations of simulation, sampling every 100th descent graph for LOD score estimation:

    swift -p east.ped -m east.map -d east.dat -o results.txt -c 4 -b 1000000 -i 1000000 -x 100

This command performs 4 separate runs and, for each run, outputs a log file starting with the prefix "log" that can be used as input to the CODA R package (see next subsection):

    swift -p east.ped -m east.map -d east.dat -o results.txt -c 4 -R 4 --trace --traceprefix log

#### Using CODA R package to perform diagnostics

This is an example to perform convergence diagnostics on the output given by the previous command in R using the CODA package. (further details about the interpretation of plots can be found on the web, for [example](http://www.johnmyleswhite.com/notebook/2010/08/29/mcmc-diagnostics-in-r-with-the-coda-package/)):

    # install package if not present
    # install.packages('coda')

    library(coda)

    chain0 <- read.table('log.ped1.run0', header=T)
    chain1 <- read.table('log.ped1.run1', header=T)
    chain2 <- read.table('log.ped1.run2', header=T)
    chain3 <- read.table('log.ped1.run3', header=T)

    chains <- mcmc.list(mcmc(chain0$likelihood), mcmc(chain1$likelihood), mcmc(chain2$likelihood), mcmc(chain3$likelihood))

    plot(chains)
    gelman.diag(chains)
    gelman.plot(chains)

# Options

    Usage: ./swift [OPTIONS] -p pedfile -m mapfile -d datfile
           ./swift [OPTIONS] -p pedfile --elod

    Input files:
      -p pedfile, --pedigree=pedfile
      -m mapfile, --map=mapfile
      -d datfile, --dat=datfile

    Output files:
      -o outfile, --output=outfile            (default = 'swiftlink.out')

    MCMC options:
      -i NUM,     --iterations=NUM            (default = 50000)
      -b NUM,     --burnin=NUM                (default = 50000)
      -s NUM,     --sequentialimputation=NUM  (default = 1000)
      -x NUM,     --scoringperiod=NUM         (default = 10)
      -l FLOAT,   --lsamplerprobability=FLOAT (default = 0.5)
      -n NUM,     --lodscores=NUM             (default = 5)
      -R NUM,     --runs=NUM                  (default = 1)

    MCMC diagnostic options:
      -T,         --trace
      -P PREFIX,  --traceprefix=PREFIX        (default = 'trace')

    ELOD options:
      -e          --elod
      -f FLOAT    --frequency=FLOAT           (default = 1.0e-04)
      -w FLOAT    --separation=FLOAT          (default = 0.0500)
      -k FLOAT,FLOAT,FLOAT --penetrance=FLOAT,FLOAT,FLOAT(default = 0.00,0.00,1.00)
      -u NUM      --replicates=NUM            (default = 1000000)

    Runtime options:
      -c NUM,     --cores=NUM                 (default = 1)
      -g,         --gpu

    Misc:
      -X,         --sexlinked
      -a,         --affectedonly
      -q NUM,     --peelseqiter=NUM           (default = 1000000)
      -r seedfile,--randomseeds=seedfile
      -v,         --verbose
      -h,         --help

