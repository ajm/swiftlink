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

 public:
	HaplotypeWriter(string fname, DescentGraph* d, Pedigree* p, GeneticMap* g) 
		: filename(fname), dg(d), ped(p), map(g) {}
	~HaplotypeWriter() {}

	bool write();
};

#endif

