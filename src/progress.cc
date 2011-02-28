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
        increments_count / increments_per_percent, 
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
        
    if((increments_count % increments_per_percent) == 0) {
        update_progress();
    }
}

void Progress::finish() { 
    complete(); 
}

