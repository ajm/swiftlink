#include <vector>
#include <sstream>

#include "trait.h"
#include "genotype.h"
#include "genetic_map.h"
#include "disease_model.h"
#include "linkage_parser.h"

// see : http://linkage.rockefeller.edu/soft/linkage/sec2.6.html
// for file format information

bool LinkageParser::unsupported(string s) {
	fprintf(stderr, "Apologises: \"%s\" not supported\n", s.c_str());
	return false;
}

bool LinkageParser::general_information(const int linenum) {
	switch(linenum) {
		case 0 :
			return set_numloci(tokens[0]) and \
				   set_sexlinked(tokens[2]) and \
				   set_programcode(tokens[3]);
		case 1 :
			return ensure_false(tokens[0], "mutation mode") and \
				   ensure_false(tokens[3], "linkage disequilibrium");
		case 2 :
			return check_num_loci();
		
		default :
			break;
	}
	
	return false;
}

bool LinkageParser::description_of_loci(const int linenum) {
	if(marker_linenum == 0) {
		if(not read_marker_code(tokens[0])) {
			return false;
		}
	}

	switch(marker_code) {
	
		case QUANTITATIVE_VARIABLE :
			return read_quantitative_variable();
			
		case AFFECTION_STATUS :
			return read_affection_status();
			
		case BINARY_FACTOR :
			return read_binary_factor();
			
		case NUMBERED_ALLELES :
			return read_numbered_allele();
			
		default :
			break;
	}
	
	// can't happen?
	fprintf(stderr, "error: %s, line %d: bad marker code %d, should be 0 - 3\n", 
				filename.c_str(), linenum, marker_code);
	return false;
}

bool LinkageParser::information_on_recombination(const int linenum) {
	vector<double>* tmp;
	
	// ignore lines 1 & 3
	// line 2 is recombination fractions
	switch(recomb_linenum) {
		case 0 : // ignore
			break;
			
		case 1 :
			tmp = get_doubles();
			if(!tmp) {
				fprintf(stderr, "error: %s, line %d: bad recombination fraction\n", 
								filename.c_str(), linenum);
				recomb_linenum++;
				return false;
			}
			
			// XXX : value zero is often for the trait marker
			// and is not at a single position in the genetic map
			// for multi-point analysis 
			for(unsigned i = 1; i < tmp->size(); ++i) {
				map.add_theta((*tmp)[i]);
			}

			delete tmp;
			break;
			
		case 2 : // ignore
			break;
			
		default :
			return false;
	}

	recomb_linenum++;

	return true;
}

bool LinkageParser::set_numloci(const string s) {
	stringstream ss(s);
	
	if((ss >> number_of_loci).fail()) {
		fprintf(stderr, "error: %s, line 1: number of loci \"%s\" is not an integer\n",
				filename.c_str(), s.c_str());
		return false;
	}
	
	return true;
}

bool LinkageParser::set_sexlinked(const string s) {
	if(s.length() == 1) {
		switch(s[0]) {
			case '0' :
				dis.set_sexlinked(false);
				return true;
			case '1' :
				dis.set_sexlinked(true);
				return true;
			default :
				break;
		}
	}

	fprintf(stderr, "error: %s, line 1: sex-linked field \"%s\" must be either "
			"0 or 1\n", filename.c_str(), s.c_str());

	return false;
}

bool LinkageParser::set_programcode(const string s) {
	stringstream ss(s);
	
	if((ss >> program_code).fail()) {
		fprintf(stderr, "error: %s, line 1: program code \"%s\" is not an integer\n",
				filename.c_str(), s.c_str());
		return false;
	}
	
	return true;
}

bool LinkageParser::ensure_false(const string s, const string msg) {
	stringstream ss(s);
	bool b;
	
	if((ss >> b).fail()) {
		fprintf(stderr, "error: %s: %s should be set to 0\n", 
			filename.c_str(), msg.c_str());
		return false;
	}

	if(b) {
		fprintf(stderr, "error: %s: %s should be set to 0", 
			filename.c_str(), msg.c_str());
	}

	return not b;
}

bool LinkageParser::check_num_loci() {
	if(int(tokens.size()) == number_of_loci)
		return true;

	fprintf(stderr, "error: %s: number of loci from line 1 (%d) and "
				"length of loci ordering from line 3 (%d) should match\n",
				filename.c_str(), number_of_loci, int(tokens.size()));
	
	return false;
}

vector<double>* LinkageParser::get_doubles() {
	vector<double>* tmp = new vector<double>();
	double tmp2;

	for(int i = 0; i < int(tokens.size()); ++i) {
		stringstream ss(tokens[i]);
		if((ss >> tmp2).fail()) {
			delete tmp;
			return NULL;
		}
		tmp->push_back(tmp2);
	}

	return tmp;
}

bool LinkageParser::read_quantitative_variable() {
	return unsupported("quantitative variable");
}

bool LinkageParser::read_binary_factor() {
	return unsupported("binary factor");
}

//
// XXX first two lines of all markers are the same
//
bool LinkageParser::read_abstract_marker() {
	vector<double>* af;
	stringstream ss;

	switch(marker_linenum) {
		case 0 :
			ss << tokens[1];
			if((ss >> marker_alleles).fail()) {
				fprintf(stderr, "error: %s, line %d: first line of marker"
							"description should have been \"%s  N\" where N is "
							"the number of alleles\n",
							filename.c_str(), linenum, tokens[0].c_str());
				return false;
			}
			marker_linenum++;
			return true;
			
		case 1 :
			af = get_doubles();

			if(marker_alleles != int(af->size())) {
				fprintf(stderr, "error: %s, line %d: expected %d alleles, but read %d\n",
							filename.c_str(), linenum, marker_alleles, int(af->size()));
			    return false;
			}
			
			// XXX assume everything is a SNP
			// or a trait for now, at least
			if(marker_alleles != 2) {
			    fprintf(stderr, "error: %s, line %d: this program is only designed to handle markers with 2 alleles\n",
							filename.c_str(), linenum);
			    return false;
			}
			
			// XXX only suitable for SNPs and traits
			marker_freq = (*af)[0];
			
			delete af;
			marker_linenum++;
			
			return true;
			
		default :
			break;
	}
	
	return false;
}

//
// XXX can only have a single liability class at the moment
//
bool LinkageParser::read_affection_status() {
	stringstream ss;
	vector<double>* af;
	bool tmp;
	
	if(total_markers_read() != 0) {
	    fprintf(stderr, "error: %s, line %d: trait marker must be the first marker in linkage file\n",
							filename.c_str(), linenum);
	    return false;
	}

	switch(marker_linenum) {
		case 0 :
		case 1 :
			return read_abstract_marker();
		
		case 2 :
			ss << tokens[0];
			if((ss >> liability_classes).fail()) {
				fprintf(stderr, "error: %s, line %d: number of liability classes"
							" must be an integer (read \"%s\")\n", 
							filename.c_str(), linenum, tokens[0].c_str());
				return false;
			}

			if(liability_classes < 1) {
				fprintf(stderr, "error: %s, line %d: must specify at least one" 
							" liability class\n",
							filename.c_str(), linenum);
				return false;
			}

			// XXX
			if(liability_classes != 1) {
				unsupported("more than one liability class");
			}

			marker_linenum++;

			return true;
		
		default :
			af = get_doubles();
			
			// store liability classes?
			if(int(af->size()) == 3) {
				for(int i = 0; i < 3; ++i)
					dis.set_penetrance((*af)[i], static_cast<enum unphased_trait>(i));
				
				tmp = true;
			}
			else {
				fprintf(stderr, "error: %s, line %d: liability classes should "
							"contain three numbers (read %d)\n",
							filename.c_str(), linenum, int(af->size()));
				tmp = false;
			}

			// seems like they are in the opposite order to numbered alleles (?)
			dis.set_freq(1 - marker_freq);
			
			delete af;
			marker_end();
			return tmp;
	}

	return false;
}

bool LinkageParser::read_numbered_allele() {
	if(! read_abstract_marker()) 
		return false;

	if(marker_linenum >= 2) {
		map[markers_read[NUMBERED_ALLELES]].set_minor_freq(marker_freq);
		marker_end();
	}

	return true;
}

bool LinkageParser::read_marker_code(const string s) {
	//printf("read_marker_code: %s\n", s.c_str());

	if(s.length() != 1) {
		fprintf(stderr, "error: %s, line %d: bad marker code %s, should be 0 - 3\n", 
				filename.c_str(), linenum, s.c_str());
		return false;
	}

	stringstream ss(s);
	if((ss >> marker_code).fail()) {
		fprintf(stderr, "error: %s, line %d: bad marker code %s, should be 0 - 3\n", 
				filename.c_str(), linenum, s.c_str());
		return false;
	}

	if((marker_code >= 0) and (marker_code < 4))
		return true;

	fprintf(stderr, "error: %s, line %d: bad marker code %s, should be 0 - 3\n", 
				filename.c_str(), linenum, s.c_str());
	return false;
}

int LinkageParser::total_markers_read() {
	int tmp = 0;
	
	for(int i = 0; i < 4; ++i) {
		tmp += markers_read[i];
    }
		
	return tmp;
}

void LinkageParser::marker_end() { 
	marker_linenum = 0; 
	markers_read[marker_code]++; 
}

bool LinkageParser::parse_line(const int linenum, const string line) {
	tokenise(line);
	
	if(linenum < 3) {
		return general_information(linenum);
	}
    
	if(total_markers_read() < number_of_loci) {
		return description_of_loci(linenum);
	}
	
	return information_on_recombination(linenum);
}

bool LinkageParser::parse_end() {
	// handle incompletely read markers
	if(marker_linenum != 0) {
	    fprintf(stderr, "error: %s: unexpected EOF\n", filename.c_str());
	    return false;
	}
	    
	if(markers_read[AFFECTION_STATUS] != 1) {
	    fprintf(stderr, "error: %s: did not read trait marker\n", 
				filename.c_str());
	    return false;
	}
	
	dis.finish_init();
	
	return true;
}

