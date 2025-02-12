#ifndef TRACK_GEN_HPP
#define TRACK_GEN_HPP

#include <string>
#include <vector>
#include <sstream>

extern "C" {
    #include "../inc/fields_nate.h"
}

std::vector<double> randomPointTrapOptimum(trace tr);

// Reads a PEN track from a delimited string.
std::vector<double> readPENTrack(const std::string &s, char delimiter);

// Generates an initial state when a PEN file is not used.
std::vector<double> TRACKGENERATOR(const trace &tr);

#endif // TRACK_GEN_HPP