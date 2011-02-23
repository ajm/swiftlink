#ifndef LKG_PROGRESS_H_
#define LKG_PROGRESS_H_

using namespace std;

#include <cstdio>
#include <string>


class Progress {
    
    string label;
    string blue_bg_colour_start;
    string black_bg_colour_start;
    string colour_end;
    unsigned int increments_total;
    unsigned int increments_count;
    unsigned int increments_per_percent;

    void update_progress() {
        printf("\r%s%s %2d%% %s", 
            blue_bg_colour_start.c_str(), 
            label.c_str(), 
            increments_count / increments_per_percent, 
            colour_end.c_str());
        fflush(stdout);
    }

    void complete() {
        printf("\r%s%s done %s\n",
            black_bg_colour_start.c_str(),
            label.c_str(),
            colour_end.c_str());
    }

 public:
    Progress(string s, unsigned int increments) 
        : label(s), increments_total(increments), increments_count(0) {
        increments_per_percent = increments_total / 100;

        blue_bg_colour_start = "\033[44m";
        black_bg_colour_start = "\033[40m\033[37m";
        colour_end = "\033[0m";
    }

    void start() {
        update_progress();
    }

    void increment() {
        increments_count++;
        
        if((increments_count % increments_per_percent) == 0) {
            update_progress();
        }
    }

    // XXX add a variadic, format string args
    void finish() { 
        complete(); 
    }
};

#endif

