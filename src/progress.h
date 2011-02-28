#ifndef LKG_PROGRESS_H_
#define LKG_PROGRESS_H_

using namespace std;

#include <cstdio>
#include <string>

#define PROGRESS_COLOUR "\033[44m"
#define ERROR_COLOUR    "\033[31m"
#define DONE_COLOUR     "\033[40m\033[37m"
#define END_COLOUR      "\033[0m"

class Progress {
    
    string label;
    unsigned int increments_total;
    unsigned int increments_count;
    unsigned int increments_per_percent;

    void update_progress();
    void complete();

 public:
    Progress(string s, unsigned int increments);

    void start();
    void increment();
    void finish();
    void error(const string& err);

};

#endif

