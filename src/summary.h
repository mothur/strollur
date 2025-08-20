//
//  summary.h
//
//  rdataset pacakge
//
//  Modified by Sarah Westcott on 5/14/25.
//  Copyright © 2025 Schloss Lab. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include "dataset.h"
#include "../inst/include/rdataset.h"

class Summary {

public:

    Summary(int p);
    ~Summary() = default;

    // report[0] = starts, report[1] = ends, report[2] = lengths,
    // report[3] = ambigs, report[4] = homopolymers, report[5] = num_ns
    Rcpp::DataFrame summarizeFasta(vector<vector<int>> report,
                                   vector<int> counts);

    Rcpp::DataFrame summarize(Rcpp::DataFrame, vector<int> counts,
                              vector<string> dfNames);
private:
    int processors;
    long long total, numUniques;

    // dataframe column names -> map of values and counts
    // "starts" -> (position -> number of seqs with that start position)
    // "ends" -> (position -> number of seqs with that end position)
    map<string, map<double, long long> > results;

    void createThreadsFasta(vector<vector<int>>& fasta,
                                          vector<int>& counts);

    void createThreadsReport(vector<vector<double>>& report,
                            vector<int>& counts, vector<string>&);

    vector<double> getValues(map<double, long long>& positions);
    vector<long long> getDefaults();
};
/******************************************************************************/

#endif /* summary_hpp */
