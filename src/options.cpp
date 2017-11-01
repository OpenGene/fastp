#include "options.h"

Options::Options(){
    in1 = "";
    in2 = "";
    out1 = "";
    out2 = "";
    thread = 1;
    compression = 3;
}

void Options::init(){
}

bool Options::isPaired(){
    return in2.length() > 0;
}