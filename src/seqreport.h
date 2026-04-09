//
//  seqreport.h
//
//  strollur package
//
//  Created by Sarah Westcott on 5/15/25.
//  Copyright © 2025 Schloss Lab. All rights reserved.
//

#ifndef SEQREPORT_H
#define SEQREPORT_H

/*
 * This class will calculate the starts, ends, lengths, ambigs, homopolyers and
 * numns for fasta sequences.
 */

#include "dataset.h"

class SeqReport {

public:

    SeqReport() { runDegap = true; }
    ~SeqReport() = default;

    // start, end, numbase, ambig, polymer, numn
    vector<int> getReport(string seq);

    // results[0] - starts, results[1] - ends, results[2] - numbases,
    // results[3] - ambigs, results[4] - polymers, results[5] - numns
    vector<vector<int>> getReport(vector<string>& seq);

    void addReports(vector<string>& seqs, vector<int>& starts, vector<int>& ends,
                   vector<int>& numbases, vector<int>& ambigs,
                   vector<int>& polymers, vector<int>& numns);

    int getStart(const string& seq) const;
    int getEnd(const string& seq) const;
    int getNumAmbigs(const string& seq) const;
    int getLongestHomopolymer(const string& seq) const;
    int getNumbases(const string& seq) const;
    int getNumns(const string& seq) const;

private:

    // used to save time when running getReport.
    bool runDegap;

    string degapSeq(const string& seq) const {

        string unaligned = seq;

        if(seq.find_first_of('.') != string::npos ||
           seq.find_first_of('-') != string::npos) {
            unaligned = "";
            for(size_t j = 0; j < seq.length(); j++) {
                if(isalpha(seq[j]))	{	unaligned += seq[j];	}
            }
        }

        return unaligned;
    }
};

#endif
