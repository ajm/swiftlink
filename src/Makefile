CXX	:= g++
CXXFLAGS := -g -Wall -Werror -Wextra -pedantic -std=c++98 -pipe 
# -Weffc++ # too hardcore...
INCLUDES := -I.
LIBS := 

OBJECTS	:= main.o \
    program.o \
    linkage_program.o \
    pedigree_parser.o \
    linkagefileparser.o \
    mapfileparser.o \
	genotype.o \
	pedigree.o \
	person.o \
	disease_model.o \
	peel_sequence_generator.o \
	rfunction.o
    

app: $(OBJECTS)
	$(CXX) -o $@ $(LIBS) $(OBJECTS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) app

