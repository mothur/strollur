//
//  summary.h
//
//  strollur pacakge
//
//  Modified by Sarah Westcott on 5/14/25.
//  Copyright © 2025 Schloss Lab. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include "dataset.h"
#include "../inst/include/strollur.h"

class Summary {

public:

    Summary(int p);
    ~Summary() = default;

    // report[0] = starts, report[1] = ends, report[2] = lengths,
    // report[3] = ambigs, report[4] = homopolymers, report[5] = num_ns
    Rcpp::DataFrame summarizeFasta(vector<vector<int>> report,
                                   vector<float> counts);

    Rcpp::DataFrame summarize(Rcpp::DataFrame, vector<float> counts,
                              vector<string> dfNames);
private:
    int processors;
    double total, numUniques;

    // dataframe column names -> map of values and counts
    // "starts" -> (position -> number of seqs with that start position)
    // "ends" -> (position -> number of seqs with that end position)
    map<string, map<double, double> > results;

    void createThreadsFasta(const vector<vector<int>>& fasta,
                                          const vector<float>& counts);

    void createThreadsReport(const vector<vector<double>>& report,
                            const vector<float>& counts, const vector<string>&);

    vector<double> getValues(map<double, double>& positions) const;
    vector<double> getDefaults() const;
};
/******************************************************************************/

#endif /* summary_hpp */
