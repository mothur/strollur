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

const float EPSILON = 0.00001f;

/**********************************************************************/
static inline bool isacgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

static bool inline isgap(char c) {
    return c == '-' || c == '.';
}

static bool inline setContains(const set<string>& t, string tag) {
    return (t.find(tag) != t.end());
}

static bool inline vectorContains(const vector<string>& t, string tag) {
    return (std::find(t.begin(), t.end(), tag) != t.end());
}

static bool inline isZero(const float& t) {
    return ((fabs(t - 0.0f)) < EPSILON);
}

static bool inline isEqual(const float& t, const float value) {
    return (isZero(t-value));
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
vector<T> setDiff(vector<T> x, vector<T> y) {

    // sort so we can use set_difference
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());

    vector<T> results;
    set_difference(x.begin(), x.end(), y.begin(), y.end(),
                   back_inserter(results));

    return results;
}
/**********************************************************************/
template<typename T>
bool addNextColumn(Rcpp::DataFrame& data, map<int, vector<T>>& y, int index) {

    auto it = y.find(index);

    if (it != y.end()) {
        data.push_back(it->second);
        return true;
    }
    return false;
}
/**********************************************************************/
// parse s by delim store in result. return numItems
template <typename Out>
int split(const string &s, char delim, Out result) {
    istringstream iss(s);
    string item;
    int numItems = 0;
    while (getline(iss, item, delim)) {
        if (!item.empty()) { //ignore white space
            *result++ = item;
            numItems++;
        }
    }
    return numItems;
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
template<typename T, typename T2>
vector<T> getValues(const map<T2, T>& x) {
    vector<T> results(x.size());

    if (x.size() == 0) { return results; }

    int index = 0;
    for (auto it = x.begin(); it != x.end(); it++) {
        results[index] = (it->second);
        index++;
    }

    return results;
}
/**********************************************************************/
template<typename T, typename T2>
vector<T> getKeys(const map<T, T2>& x) {
    vector<T> results(x.size());

    if (x.size() == 0) { return results; }

    int index = 0;
    for (auto it = x.begin(); it != x.end(); it++) {
        results[index] = (it->first);
        index++;
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
T sum(const vector<T>& x) {
    return accumulate(x.begin(), x.end(), 0);
}
/**********************************************************************/
template<typename T>
bool allIdentical(const vector<T>& x, T& value) {

    if (!x.empty()) {
        value = x[0];
        for (int i = 0; i < x.size(); i++) {
            if (x[i] != value) {
                return false;
            }
        }

    }
    return true;
}
/**********************************************************************/
static inline bool allBlank(const vector<string>& x) {
    string value = "";
    if (allIdentical(x, value)) {
        if (value == "") {
            return true;
        }
    }
    return false;
}
/**********************************************************************/
template<typename T>
bool identical(vector<T> x, vector<T> y) {

    sort(x.begin(), x.end());
    sort(y.begin(), y.end());

    return (x == y);
}
/**********************************************************************/
template<typename T>
void sum(vector<T>& x, vector<T>& y) {
    transform (x.begin(), x.end(),
               y.begin(), x.begin(), std::plus<T>());
}
/**********************************************************************/
template<typename T>
void subtract(vector<T>& x, vector<T>& y) {
    transform (x.begin(), x.end(),
               y.begin(), x.begin(), std::minus<T>());
}
/**********************************************************************/
template<typename T>
vector<T> unique(const vector<T>& x) {

    vector<T> uniqueX;

    if (x.size() == 0) { return uniqueX; }

    set<T> s;
    unsigned size = x.size();
    for( unsigned i = 0; i < size; ++i ) {
        s.insert( x[i] );
    }

    uniqueX.assign( s.begin(), s.end() );

    return uniqueX;
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
static bool inline isTrue(vector<bool> c) {

    if (c.size() == 1) { return c[0]; }

    set<bool> x = toSet(c);
    if (x.size() > 1) {
        return false;
    }else if (x.size() == 1) {
        return (*x.begin() == true);
    }
    return true;
}
/**********************************************************************/
static bool inline isFalse(vector<bool> c) {

    if (c.size() == 1) { return (c[0] == false); }

    set<bool> x = toSet(c);

    if (x.size() > 1) {
        return false;
    }else if (x.size() == 1) {
        return (*x.begin() == false);
    }
    return true;
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
class BadConversion : public runtime_error {
public:
    BadConversion(const string& s) : runtime_error(s){ }
};
/**********************************************************************/
template<typename T>
void convert(const string& s, T& x, bool failIfLeftoverChars = true){

    istringstream i(s);
    char c;
    if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
        throw BadConversion(s);

}
/**********************************************************************/

#endif
