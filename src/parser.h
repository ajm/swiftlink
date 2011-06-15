#ifndef LKG_PARSER_H_
#define LKG_PARSER_H_

using namespace std;

#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <istream>
#include <fstream>


class Parser {

 protected:
	string filename;
	bool die_on_error;
	int linenum;
	ifstream f;
	vector<string> tokens;

	void tokenise(string s) {
		stringstream ss(s);
        
        istream_iterator<string> i(ss);
        istream_iterator<string> end;
        
        this->tokens = vector<string>(i,end);
	}

 public:
	Parser(const string fn, bool die_on_err) 
		: filename(fn), die_on_error(die_on_err), linenum(0), f(), tokens() {}
	virtual ~Parser() {}
	virtual bool parse_line(const int linenum, const string s) = 0;
	virtual bool parse_end() { return true; }
	bool parse() {
		string line;
		bool noerror = true;
        string::size_type index;    // OSX bitches about signed vs. unsigned 
									// comparisons and if they are both unsigned
									// says this is always true due to limied 
									// range of data type  :-S

		f.open(filename.c_str());

		if(!f) {
			fprintf(stderr, "error: file not found: %s\n", filename.c_str());
			return false;
		}

		while(! getline(f, line).eof()) {
			// terminate string from comment char
			index = line.find("#");
			if(index != string::npos) {
				line.erase(line.begin() + index, line.end());
			}

			// forget about empty lines
			if(line.length() == 0)
				continue;

			if(! parse_line(linenum, line)) {
				noerror = false;
			}

			if(die_on_error and not noerror)
				break;

			linenum++;
		}

		f.close();

		if(! parse_end()) {
			noerror = false;
		}
		
		return noerror;
	}
};

// simple example of how you would use parser 
// abstract class
class EchoParser : public Parser {

 public :
	EchoParser(const string fn) 
		: Parser(fn, false) {}

	bool parse_line(const int linenum, const string s) {
		tokenise(s);
		
        printf("line %d\n", linenum);
		for(int i = 0; i < int(tokens.size()); ++i) {
			printf("%s:", tokens[i].c_str());
		}
        printf("\n");
		
		return true;
	}
};

#endif

