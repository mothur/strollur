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

const vector<int> nullIntVector;

/**********************************************************************/
struct pieceOfWork {
    double start;
    double end;
    pieceOfWork(double i, double j) : start(i), end(j) {}
    pieceOfWork() { start = 0; end = 0; }
    ~pieceOfWork() {}
};
/**********************************************************************/
vector<pieceOfWork> divideWork(double numItems, int& numProcessors) {
    // divide work between processors
    vector<pieceOfWork> work;

    if (numItems < numProcessors) { numProcessors = numItems; }
    size_t startIndex = 0;

    for (size_t remainingProcessors = numProcessors; remainingProcessors > 0;
    remainingProcessors--) {

        //case for last processor
        size_t numToProcess = numItems;
        if (remainingProcessors != 1) {
            numToProcess = ceil(numItems / remainingProcessors);
        }
        work.push_back(pieceOfWork(startIndex, (startIndex+numToProcess)));
        startIndex += numToProcess;
        numItems -= numToProcess;
    }
    return work;
}
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
