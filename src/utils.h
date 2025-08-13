#ifndef utils_h
#define utils_h

#include "dataset.h"

class Utils {

public:

    Utils() = default;
    ~Utils() = default;

    // taxonomy helpers
    //void addUnclassifieds(string& taxon, int maxlevel);
    void addUnclassifieds(vector<string>& taxons,
                          vector<int>& confidences, int maxlevel);
    bool findTaxon(vector<string> tax, vector<string> stax);
    int removeConfidence(string& taxon);
    vector<int> removeConfidences(vector<string>& taxons);
    bool searchTax(vector<string> userTaxons,
                   vector<int> userConfidences,
                   vector<bool> taxonsHasConfidence,
                   vector< vector<string> > searchTaxons,
                   vector< vector<int> > searchConfidenceThresholds);

    // paralell processing helper
    vector<pieceOfWork> divideWork(double numItems, int& numProcessors);

    // numeric helpers
    bool isPositiveNumeric(string s);

private:

    string removeQuotes(string tax);
    void removeQuotes(vector<string>& tax);

};

#endif /* utils_h */

