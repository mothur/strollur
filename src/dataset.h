#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

// containers
#include <vector>
#include <map>
#include <string>
#include <string.h>

// Rcpp
#include <Rcpp.h>
#include <RcppThread.h>
#include <cli/progress.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

using namespace std;

/**********************************************************************/
static inline bool isacgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

static bool inline isgap(char c) {
    return c == '-' || c == '.';
}
/**********************************************************************/
struct pieceOfWork {
    double start;
    double end;
    pieceOfWork(double i, double j) : start(i), end(j) {}
    pieceOfWork() { start = 0; end = 0; }
    ~pieceOfWork() {}
};
/**********************************************************************/
template<typename T>
set<T> toSet(const vector<T>& x) {
    set<T> results;

    if (x.size() == 0) { return results; }

    for (int i = 0; i < x.size(); i++ ) {
        results.insert(x[i]);
    }

    return results;
}
/**********************************************************************/
template<typename T>
vector<T> select(const vector<T>& x, const vector<bool>& filter) {
    vector<T> results;

    if (x.size() == 0) { return results; }

    for (int i = 0; i < x.size(); i++ ) {
        if (filter[i]) {
            results.push_back(x[i]);
        }
    }

    return results;
}
/**********************************************************************/
template<typename T>
vector<T> toVector(const set<T>& x) {
    vector<T> results;

    if (x.size() == 0) { return results; }

    for (auto it = x.begin(); it != x.end(); it++) {
        results.push_back(*it);
    }

    return results;
}
/**********************************************************************/
inline string toString(const set<string>& x, char delim) {
    string result = "";

    if (x.size() == 0) { return result; }

    for (auto it = x.begin(); it != x.end(); it++) {
        result += delim + *it;
    }
    result = result.substr(1);

    return result;
}
/**********************************************************************/
inline string toString(const bool& x) {
    if (x) {
        return "TRUE";
    }else{
        return "FALSE";
    }
}
/**********************************************************************/
inline string toString(const set<string>& x, string delim) {
    string result = "";

    if (x.size() == 0) { return result; }

    for (auto it = x.begin(); it != x.end(); it++) {
        result += delim + *it;
    }
    result = result.substr(delim.length());

    return result;
}
/**********************************************************************/
template<typename T>
string toString(const T&x) {
    stringstream output;
    output << x;
    return output.str();
}
/**********************************************************************/
template<typename T>
string toString(const T&x, int i) {
    stringstream output;

    output.precision(i);
    output << std::fixed << x;

    return output.str();
}
/**********************************************************************/
template<typename T>
string toString(const vector<T>& x, char delim) {
    string result = "";

    if (x.size() == 0) { return result; }

    result = toString(x[0]);

    for (int i = 1; i < x.size(); i++) {
        result += delim + toString(x[i]);
    }

    return result;
}
/**********************************************************************/

#endif
