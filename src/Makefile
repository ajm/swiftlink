CXX	:= g++
CXXFLAGS := -g -Wall -Werror -Wextra -pedantic -std=c++98 -pipe 
# -Weffc++ # too hardcore...
INCLUDES := -I.
LIBS := 

OBJECTS	:= test_peeling.o \
	pedigree.o \
	person.o \
	pedigree_parser.o \
	disease_model.o \
	genotype.o \
	peel_sequence_generator.o \
	rfunction.o

test_peeling: $(OBJECTS)
	$(CXX) -o $@ $(LIBS) $(OBJECTS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) test_peeling

