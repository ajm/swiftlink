#ifndef LKG_HAPLOTYPEWRITER_H_
#define LKG_HAPLOTYPEWRITER_H_

using namespace std;

#include <string>

class DescentGraph;
class Pedigree;
class GeneticMap;

class HaplotypeWriter {
	string filename;
	DescentGraph* dg;
	Pedigree* ped;
	GeneticMap* map;
	bool verbose;

 public:
	HaplotypeWriter(DescentGraph* d, Pedigree* p, GeneticMap* g, const char* fname, bool verbose) 
		: filename(fname), dg(d), ped(p), map(g), verbose(verbose) {}
	~HaplotypeWriter() {}

	bool write();
};

#endif

