#ifndef LKG_LINKAGEWRITER_H_
#define LKG_LINKAGEWRITER_H_

using namespace std;

#include <string>

class Peeler;
class GeneticMap;

class LinkageWriter {
	GeneticMap* map;
	Peeler* peeler;
	string filename;
	bool verbose;

 public:
	LinkageWriter(GeneticMap* g, Peeler* peel, const char* filename, bool verbose) 
		: map(g), peeler(peel), filename(filename), verbose(verbose) {}
	~LinkageWriter() {}

	bool write();
};

#endif

