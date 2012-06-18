#ifndef LKG_LINKAGEWRITER_H_
#define LKG_LINKAGEWRITER_H_

using namespace std;

#include <string>
#include <vector>

class Peeler;
class GeneticMap;
class LODscores;

class LinkageWriter {
	GeneticMap* map;
	string filename;
	bool verbose;

 public:
	LinkageWriter(GeneticMap* g, string filename, bool verbose) : 
	    map(g), 
	    filename(filename), 
	    verbose(verbose) {}
	
	~LinkageWriter() {}

    LinkageWriter(const LinkageWriter& rhs) :
        map(rhs.map),
        filename(rhs.filename),
        verbose(rhs.verbose) {}
        
    LinkageWriter& operator=(const LinkageWriter& rhs) {
        if(&rhs != this) {
            map = rhs.map;
            filename = rhs.filename;
            verbose = rhs.verbose;
        }
        
        return *this;
    }

	bool write(vector<LODscores*>& all_scores);
};

#endif

