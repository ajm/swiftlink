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
	LinkageWriter(GeneticMap* g, Peeler* peel, const char* filename, bool verbose) : 
	    map(g), 
	    peeler(peel), 
	    filename(filename), 
	    verbose(verbose) {}
	
	~LinkageWriter() {}

    LinkageWriter(const LinkageWriter& lw) :
        map(lw.map),
        peeler(lw.peeler),
        filename(lw.filename),
        verbose(lw.verbose) {}
        
    LinkageWriter& operator=(const LinkageWriter& rhs) {
        if(&rhs != this) {
            map = rhs.map;
            peeler = rhs.peeler;
            filename = rhs.filename;
            verbose = rhs.verbose;
        }
        
        return *this;
    }

	bool write();
};

#endif

