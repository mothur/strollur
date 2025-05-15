//
//  summary.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/17.
//  Copyright © 2024 Schloss Lab. All rights reserved.
//

#include "summary.h"
#include "utils.h"

//**************************************************************************
Summary::Summary(int p) {
    processors = p;
    total = 0;
    numUniques = 0;
    hasCount = false;
}
//**************************************************************************
vector<long long> Summary::getDefaults() {
        vector<long long> locations;

        //number of sequences at 2.5%
        long long ptile0_25	= 1+(long long)(total * 0.025);
        //number of sequences at 25%
        long long ptile25		= 1+(long long)(total * 0.250);
        long long ptile50		= 1+(long long)(total * 0.500);
        long long ptile75		= 1+(long long)(total * 0.750);
        long long ptile97_5	= 1+(long long)(total * 0.975);
        long long ptile100	= (long long)(total);

        locations.push_back(1);
        locations.push_back(ptile0_25);
        locations.push_back(ptile25);
        locations.push_back(ptile50);
        locations.push_back(ptile75);
        locations.push_back(ptile97_5);
        locations.push_back(ptile100);

        return locations;
}
//**************************************************************************
long long Summary::getValue(map<int, long long>& spots, double value) {

    long long percentage = 1+(long long)(total * value * 0.01);
    long long result = 0;
    long long totalSoFar = 0;
    long long lastValue = 0;

    //minimum
    if ((spots.begin())->first == -1) { result = 0; }
    else {result = (spots.begin())->first; }

    for (auto it = spots.begin(); it != spots.end(); it++) {
        long long value = it->first; if (value == -1) { value = 0; }
        totalSoFar += it->second;

        if (((totalSoFar <= percentage) && (totalSoFar > 1)) ||
            ((lastValue < percentage) && (totalSoFar > percentage))){
            result = value;
        }
        lastValue = totalSoFar;
    }

    return result;

}
//**************************************************************************
vector<int> Summary::getValues(map<int, long long>& positions) {
        vector<long long> defaults = getDefaults();
        vector<int> results; results.resize(7,0);
        long long meanPosition; meanPosition = 0;
        long long totalSoFar = 0;
        int lastValue = 0;

        // minimum
        if ((positions.begin())->first != -1) {
            results[0] = (positions.begin())->first;
        }

        results[1] = results[0]; results[2] = results[0];
        results[3] = results[0]; results[4] = results[0];
        results[5] = results[0];

        for (auto it = positions.begin(); it != positions.end(); it++) {

            int value = it->first; if (value == -1) { value = 0; }
            meanPosition += (value*it->second);
            totalSoFar += it->second;

            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) ||
                ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){
                results[1] = value;
            } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||
                ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) {
                results[2] = value;
            } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||
                ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {
                results[3] = value;
            } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||
                ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {
                results[4] = value;
            } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||
                ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {
                results[5] = value;
            } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {
                results[6] = value;
            } //save value
            lastValue = totalSoFar;
        }
        // maximum
        results[6] = (positions.rbegin())->first;


        double meansPosition = meanPosition / (double) total;

        // mean
        results.push_back(meansPosition);

        return results;
}

//**************************************************************************
vector<float> Summary::getValues(map<float, long long>& positions) {
        vector<long long> defaults = getDefaults();
        vector<float> results; results.resize(7,0);
        long long meanPosition; meanPosition = 0;
        long long totalSoFar = 0;
        int lastValue = 0;

        //minimum
        if ((positions.begin())->first != -1) {
            results[0] = (positions.begin())->first;
        }

        results[1] = results[0]; results[2] = results[0];
        results[3] = results[0]; results[4] = results[0];
        results[5] = results[0];

        for (auto it = positions.begin(); it != positions.end(); it++) {

            float value = it->first; if (value == -1) { value = 0; }
            meanPosition += (value*it->second);
            totalSoFar += it->second;

            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) ||
                ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){
                results[1] = value;
            } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||
                ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) {
                results[2] = value;
            } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||
                ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {
                results[3] = value;
            } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||
                ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {
                results[4] = value;
            } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||
                ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {
                results[5] = value;
            } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {
                results[6] = value;
            } //save value
            lastValue = totalSoFar;
        }
        // maximum
        results[6] = (positions.rbegin())->first;

        double meansPosition = meanPosition / (double) total;

        // mean
        results.push_back(meansPosition);

        return results;
}
//**************************************************************************
struct seqSumData {
    // fasta
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;

    //contigs
    map<int, long long> ostartPosition;
    map<int, long long> oendPosition;
    map<int, long long> oseqLength;
    map<int, long long> misMatches;
    map<int, long long> numNs;

    //align
    map<float, long long> sims;
    map<float, long long> scores;
    map<int, long long> inserts;

    vector<int> counts; // abundances for each read, all ones if no count info
    vector<vector<int>> contigs;
    vector<vector<float>> align;
    vector<int> alignInserts;
    vector<vector<int>> summary;

    double start;
    double end;

    seqSumData(){
        start = 0;
        end = 0;
    }
    // fasta and contigs summary constructor
    seqSumData(vector<vector<int>>& f, double st, double en, vector<int> nam,
               string method) {
        if (method == "fasta") {
            summary = f;
        }else if (method == "contigs") {
            contigs = f;
        }

        start = st;
        end = en;
        counts = nam;
    }

    // align summary constructor
    seqSumData(vector<vector<float>>& f, vector<int> i, double st, double en,
               vector<int> nam) {
        align= f;
        alignInserts = i;
        start = st;
        end = en;
        counts = nam;
    }
};
//**************************************************************************
Rcpp::DataFrame Summary::summarizeFasta(vector<vector<int>>& report,
                                  vector<int> counts) {
        numUniques = counts.size();
        total = numUniques;
        if (counts.size() != 0) {
            int initial_sum = 0;
            total = accumulate(counts.begin(), counts.end(), initial_sum);
        }else {
            counts.resize(numUniques, 1);
        }

        createThreadsFasta(report, counts);

         // can't pass long long, convert to string
         vector<long long> numSeqs = getDefaults();
         vector<string> wrappedNumSeqs(numSeqs.size()+1, "");
         for (int i = 0; i < numSeqs.size(); i++) {
             wrappedNumSeqs[i] = toString(numSeqs[i]);
         }

         Rcpp::DataFrame df = Rcpp::DataFrame::create(
             Rcpp::Named("starts") = getStart(),
             Rcpp::_["ends"] = getEnd(),
             Rcpp::_["nbases"] = getLength(),
             Rcpp::_["ambigs"] = getAmbig(),
             Rcpp::_["polymers"] = getHomop(),
             Rcpp::_["numns"] = getNumNs(),
             Rcpp::_["numseqs"] = wrappedNumSeqs);
        return (df);
}
//**************************************************************************

void driverSummarize(seqSumData* params) {
    // calc lengths
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisSeqLength = params->summary[0][i];
        auto it = params->seqLength.find(thisSeqLength);
        if (it == params->seqLength.end()) {
            params->seqLength[thisSeqLength] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // starts
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisStartPosition = params->summary[1][i];
        auto it = params->startPosition.find(thisStartPosition);
        if (it == params->startPosition.end()) {
            params->startPosition[thisStartPosition] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // ends
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisEndPosition = params->summary[2][i];
        auto it = params->endPosition.find(thisEndPosition);
        if (it == params->endPosition.end()) {
            params->endPosition[thisEndPosition] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // ambigs
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisAmbig = params->summary[3][i];
        auto it = params->ambigBases.find(thisAmbig);
        if (it == params->ambigBases.end()) {
            params->ambigBases[thisAmbig] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    // polymers
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisHomoP = params->summary[4][i];
        auto it = params->longHomoPolymer.find(thisHomoP);
        if (it == params->longHomoPolymer.end()) {
            params->longHomoPolymer[thisHomoP] = num; }
        else { it->second += num; } //add counts

    }

    RcppThread::checkUserInterrupt();

    // numNs
    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisNumNs = params->summary[5][i];
        auto it = params->numNs.find(thisNumNs);
        if (it == params->numNs.end()) {
            params->numNs[thisNumNs] = num; }
        else { it->second += num; } //add counts
    }
}
//**************************************************************************
void Summary::createThreadsFasta(vector<vector<int>>& summary,
                                 vector<int>& counts) {
         //divide reads between processors
         Utils util;
        vector<pieceOfWork> startEndIndexes = util.divideWork(counts.size(),
                                                              processors);

        vector<RcppThread::Thread*> workerThreads;
        vector<seqSumData*> data;

        //Launch worker threads
        for (int i = 0; i < processors-1; i++) {

            seqSumData* dataBundle = new seqSumData(summary,
                                                    startEndIndexes[i+1].start,
                                                    startEndIndexes[i+1].end,
                                                    counts, "fasta");
            data.push_back(dataBundle);

            workerThreads.push_back(new RcppThread::Thread(driverSummarize,
                                                           dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(summary, startEndIndexes[0].start,
                                                startEndIndexes[0].end,
                                                counts, "fasta");
        driverSummarize(dataBundle);
        startPosition = dataBundle->startPosition;
        endPosition = dataBundle->endPosition;
        seqLength = dataBundle->seqLength;
        ambigBases = dataBundle->ambigBases;
        longHomoPolymer = dataBundle->longHomoPolymer;
        numNs = dataBundle->numNs;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();

            for (auto it = data[i]->startPosition.begin();
                 it != data[i]->startPosition.end(); it++) {
                auto itMain = startPosition.find(it->first);
                if (itMain == startPosition.end()) { //newValue
                    startPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->endPosition.begin();
                 it != data[i]->endPosition.end(); it++) {
                auto itMain = endPosition.find(it->first);
                if (itMain == endPosition.end()) { //newValue
                    endPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->seqLength.begin();
                 it != data[i]->seqLength.end(); it++)		{
                auto itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->ambigBases.begin();
                 it != data[i]->ambigBases.end(); it++)		{
                auto itMain = ambigBases.find(it->first);
                if (itMain == ambigBases.end()) { //newValue
                    ambigBases[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->longHomoPolymer.begin();
                 it != data[i]->longHomoPolymer.end(); it++)		{
                auto itMain = longHomoPolymer.find(it->first);
                if (itMain == longHomoPolymer.end()) { //newValue
                    longHomoPolymer[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            for (auto it = data[i]->numNs.begin();
                 it != data[i]->numNs.end(); it++)		{
                auto itMain = numNs.find(it->first);
                if (itMain == numNs.end()) { //newValue
                    numNs[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }
}
//**************************************************************************
Rcpp::DataFrame Summary::summarizeContigs(vector<vector<int>>& report,
                                    vector<int> counts) {

        numUniques = counts.size();
        int initial_sum = 0;
        total = accumulate(counts.begin(), counts.end(), initial_sum);

        createThreadsContigs(report, counts);

        // can't pass long long, convert to string
        vector<long long> numSeqs = getDefaults();
        vector<string> wrappedNumSeqs(numSeqs.size()+1, "");
        for (int i = 0; i < numSeqs.size(); i++) {
            wrappedNumSeqs[i] = toString(numSeqs[i]);
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("ostarts") = getOStart(),
            Rcpp::_["oends"] = getOEnd(),
            Rcpp::_["lengths"] = getLength(),
            Rcpp::_["olengths"] = getOLength(),
            Rcpp::_["numns"] = getNumNs(),
            Rcpp::_["mismatches"] = getMisMatches(),
            Rcpp::_["numseqs"] = wrappedNumSeqs);
        return (df);
}
//**************************************************************************

void driverContigs(seqSumData* params) {
    // names = c(), lengths = c(), olengths = c(),
    //     ostarts = c(), oends = c(), mismatches = c(),
    //     numns = c(), ee = c()
    // length overlap_length overlap_start overlap_end mismatches num_ns
    // 1.       2              3               4           5       6

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisSeqLength = params->contigs[0][i];
        auto it = params->seqLength.find(thisSeqLength);
        if (it == params->seqLength.end()) {
            params->seqLength[thisSeqLength] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisOLength = params->contigs[1][i];
        auto it = params->oseqLength.find(thisOLength);
        if (it == params->oseqLength.end()) {
            params->oseqLength[thisOLength] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisOStartPosition = params->contigs[2][i];
        auto it = params->ostartPosition.find(thisOStartPosition);
        if (it == params->ostartPosition.end()) {
            params->ostartPosition[thisOStartPosition] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisOEndPosition = params->contigs[3][i];
        auto it = params->oendPosition.find(thisOEndPosition);
        if (it == params->oendPosition.end()) {
            params->oendPosition[thisOEndPosition] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisMismatches = params->contigs[4][i];
        auto it = params->misMatches.find(thisMismatches);
        if (it == params->misMatches.end()) {
            params->misMatches[thisMismatches] = num; }
        else { it->second += num; } //add counts

    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisNumNs = params->contigs[5][i];
        auto it = params->numNs.find(thisNumNs);
        if (it == params->numNs.end()) {
            params->numNs[thisNumNs] = num; }
        else { it->second += num; } //add counts
    }
}
//**************************************************************************
void Summary::createThreadsContigs(vector<vector<int>>& contigs,
                                   vector<int>& counts) {
         //divide reads between processors
         Utils util;
        vector<pieceOfWork> startEndIndexes = util.divideWork(counts.size(),
                                                              processors);

        vector<RcppThread::Thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            seqSumData* dataBundle = new seqSumData(contigs,
                                                    startEndIndexes[i+1].start,
                                                    startEndIndexes[i+1].end,
                                                    counts, "contigs");
            data.push_back(dataBundle);

            workerThreads.push_back(new RcppThread::Thread(driverContigs,
                                                           dataBundle));
        }
        seqSumData* dataBundle = new seqSumData(contigs,
                                                startEndIndexes[0].start,
                                                startEndIndexes[0].end,
                                                counts, "contigs");

        driverContigs(dataBundle);
        ostartPosition = dataBundle->ostartPosition;
        oendPosition = dataBundle->oendPosition;
        seqLength = dataBundle->seqLength;
        oseqLength = dataBundle->oseqLength;
        misMatches = dataBundle->misMatches;
        numNs = dataBundle->numNs;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();

            for (auto it = data[i]->ostartPosition.begin();
                 it != data[i]->ostartPosition.end(); it++) {
                auto itMain = ostartPosition.find(it->first);
                if (itMain == ostartPosition.end()) { //newValue
                    ostartPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->oendPosition.begin();
                 it != data[i]->oendPosition.end(); it++) {
                auto itMain = oendPosition.find(it->first);
                if (itMain == oendPosition.end()) { //newValue
                    oendPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->seqLength.begin();
                 it != data[i]->seqLength.end(); it++)		{
                auto itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->oseqLength.begin();
                 it != data[i]->oseqLength.end(); it++)		{
                auto itMain = oseqLength.find(it->first);
                if (itMain == oseqLength.end()) { //newValue
                    oseqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->misMatches.begin();
                 it != data[i]->misMatches.end(); it++)		{
                auto itMain = misMatches.find(it->first);
                if (itMain == misMatches.end()) { //newValue
                    misMatches[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->numNs.begin();
                 it != data[i]->numNs.end(); it++)		{
                auto itMain = numNs.find(it->first);
                if (itMain == numNs.end()) { //newValue
                    numNs[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }
}
//**************************************************************************
Rcpp::DataFrame Summary::summarizeAlign(vector<vector<float>>& align,
                                  vector<int>& inserts,
                                  vector<int> counts) {
        numUniques = counts.size();
        int initial_sum = 0;
        total = accumulate(counts.begin(), counts.end(), initial_sum);

        createThreadsAlign(align, inserts, counts);

        // can't pass long long, convert to string
        vector<long long> numSeqs = getDefaults();
        vector<string> wrappedNumSeqs(numSeqs.size()+1, "");
        for (int i = 0; i < numSeqs.size(); i++) {
            wrappedNumSeqs[i] = toString(numSeqs[i]);
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("search_scores") = getScores(),
            Rcpp::_["sim_scores"] = getSims(),
            Rcpp::_["longest_inserts"] = getNumInserts(),
            Rcpp::_["numseqs"] = wrappedNumSeqs);
        return (df);
}
//**************************************************************************

void driverAlign(seqSumData* params) {
    // search_scores = 0 sim_scores = 1
    // scores, sims, inserts

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        float thisSearchScore = params->align[0][i];
        auto it = params->scores.find(thisSearchScore);
        if (it == params->scores.end()) {
            params->scores[thisSearchScore] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        float thisSimScore = params->align[1][i];
        auto it = params->sims.find(thisSimScore);
        if (it == params->sims.end()) {
            params->sims[thisSimScore] = num; }
        else { it->second += num; } //add counts
    }

    RcppThread::checkUserInterrupt();

    for (int i = params->start; i < params->end; i++) {

        int num = params->counts[i];

        int thisInserts = params->alignInserts[i];
        auto it = params->inserts.find(thisInserts);
        if (it == params->inserts.end()) {
            params->inserts[thisInserts] = num; }
        else { it->second += num; } //add counts
    }
}
//**************************************************************************
void Summary::createThreadsAlign(vector<vector<float>>& align,
                                 vector<int>& alignInserts,
                                   vector<int>& counts) {

        //divide reads between processors
        Utils util;
        vector<pieceOfWork> startEndIndexes = util.divideWork(counts.size(),
                                                              processors);

        vector<RcppThread::Thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            seqSumData* dataBundle = new seqSumData(align, alignInserts,
                                                    startEndIndexes[i+1].start,
                                                    startEndIndexes[i+1].end,
                                                    counts);
            data.push_back(dataBundle);

            workerThreads.push_back(new RcppThread::Thread(driverAlign,
                                                           dataBundle));
        }
        seqSumData* dataBundle = new seqSumData(align, alignInserts,
                                                startEndIndexes[0].start,
                                                startEndIndexes[0].end,
                                                counts);

        driverAlign(dataBundle);
        sims = dataBundle->sims;
        scores = dataBundle->scores;
        inserts = dataBundle->inserts;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();

            for (auto it = data[i]->sims.begin();
                 it != data[i]->sims.end(); it++) {
                auto itMain = sims.find(it->first);
                if (itMain == sims.end()) { //newValue
                    sims[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->scores.begin();
                 it != data[i]->scores.end(); it++) {
                auto itMain = scores.find(it->first);
                if (itMain == scores.end()) { //newValue
                    scores[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (auto it = data[i]->inserts.begin();
                 it != data[i]->inserts.end(); it++)		{
                auto itMain = inserts.find(it->first);
                if (itMain == inserts.end()) { //newValue
                    inserts[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }
}
//**************************************************************************
