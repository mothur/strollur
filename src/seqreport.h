//
//  seqreport.h
//
//  rdataset package
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

    int getStart(string seq);
    int getEnd(string seq);
    int getNumAmbigs(string seq);
    int getLongestHomopolymer(string seq);
    int getNumbases(string seq);
    int getNumns(string seq);

private:

    // used to save time when running getReport.
    bool runDegap;

    string degapSeq(string seq) {

        string unaligned = seq;

        if(seq.find_first_of('.') != string::npos ||
           seq.find_first_of('-') != string::npos) {
            unaligned = "";
            for(int j = 0; j < seq.length(); j++) {
                if(isalpha(seq[j]))	{	unaligned += seq[j];	}
            }
        }

        return unaligned;
    }
};

#endif
