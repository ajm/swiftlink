using namespace std;

#include <cstdio>
#include <string>

#include "progress.h"


Progress::Progress(string s, unsigned int increments) 
    : label(s), increments_total(increments), increments_count(0) {

    increments_per_percent = increments_total / 100;
}

void Progress::update_progress() {
    printf("\r%s%s %2d%% %s", 
        PROGRESS_COLOUR, 
        label.c_str(), 
        int((increments_count / double(increments_total)) * 100), 
        END_COLOUR);
    fflush(stdout);
}

void Progress::complete() {
    printf("\r%s%s done %s\n",
        DONE_COLOUR,
        label.c_str(),
        END_COLOUR);
}

void Progress::start() {
    update_progress();
}

void Progress::increment() {

    increments_count++;
    
    if((increments_per_percent == 0) or \
        ((increments_count % increments_per_percent) == 0)) {

        update_progress();
    }
}

void Progress::finish() { 
    complete(); 
}

void Progress::error(const string& err) {
    printf("\r%s%s %s %s\n",
        ERROR_COLOUR,
        label.c_str(),
        err.c_str(),
        END_COLOUR);
}

