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

    // start, end, numbases, ambigs, polymers, numns
    vector<int> getReport(string seq) {
        vector<int> results;

        results.push_back(getStart(seq));
        results.push_back(getEnd(seq));

        seq = degapSeq(seq);
        runDegap = false;

        results.push_back(seq.length());
        results.push_back(getNumAmbigs(seq));
        results.push_back(getLongestHomopolymer(seq));
        results.push_back(getNumns(seq));

        runDegap = true;

        return results;
    }

    int getStart(string seq) {
        int startPos = 1;
        for(int i = 0; i < seq.length(); i++) {
            if((seq[i] != '.') && (seq[i] != '-')){
                startPos = i + 1;
                break;
            }
        }
        return startPos;
    }

    int getEnd(string seq) {
        int endPos = 1;
        for(int i = seq.length()-1; i >= 0; i--){
            if((seq[i] != '.') && (seq[i] != '-')){
                endPos = i + 1;
                break;
            }
        }
        return endPos;
    }

    int getNumAmbigs(string seq) {
        int ambigBases = 0;
        for(int i = 0; i < seq.length(); i++) {
            if(!isgap(seq[i]) && !isacgt(seq[i])){
                ambigBases++;
            }
        }
        return ambigBases;
    }

    int getLongestHomopolymer(string seq) {
        int longHomoPolymer = 1;
        int homoPolymer = 1;

        string unaligned = seq;

        if (runDegap) {
            unaligned = degapSeq(seq);
        }

        for(int j = 1; j < unaligned.length(); j++){
            if(unaligned[j] == unaligned[j-1]){
                homoPolymer++;
            }else{
                if(homoPolymer > longHomoPolymer){
                    longHomoPolymer = homoPolymer;
                }
                homoPolymer = 1;
            }
        }

        if(homoPolymer > longHomoPolymer){
            longHomoPolymer = homoPolymer;
        }

        return longHomoPolymer;
    }

    int getNumbases(string seq) {
        int numBases = 0;
        for(int i = 0; i < seq.length(); i++) {
            if(!isgap(seq[i])){
                numBases++;
            }
        }
        return numBases;
    }

    int getNumns(string seq) {
        int numNs = 0;
        for (int i = 0; i < seq.length(); i++) {
            if(seq[i] == 'N') { numNs++; }
        }
        return numNs;
    }

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
