using namespace std;

#include <cstdio>
#include <string>
#include <cstdarg>

#include "progress.h"


Progress::Progress(string s, unsigned int increments) : 
    label(s), 
    increments_total(increments), 
    increments_count(0),
    increments_per_percent(increments / 1000) {

    update_progress();    
}

void Progress::update_progress() {
    /*
    printf("\r%s%s %.1f%% %s", 
        PROGRESS_COLOUR, 
        label.c_str(), 
        (increments_count / double(increments_total)) * 100, 
        END_COLOUR);
    */
    printf("\r%s %.1f%%", 
        label.c_str(), 
        (increments_count / double(increments_total)) * 100);
    fflush(stdout);
}

void Progress::finish() {
    /*
    printf("\r%s%s done   %s\n",
        DONE_COLOUR,
        label.c_str(),
        END_COLOUR);
    */
    printf("\r%s done   \n",
        label.c_str());
}

void Progress::finish_msg(const char* fmt, ...) {
    printf("\r%s ", label.c_str());
    
    va_list argptr;
    va_start(argptr, fmt);
    vprintf(fmt, argptr);
    va_end(argptr);
    
    printf("   \n");
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

void Progress::error(const string& err) {
    printf("\r%s%s %s %s\n",
        ERROR_COLOUR,
        label.c_str(),
        err.c_str(),
        END_COLOUR);
}

