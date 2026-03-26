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
}
//**************************************************************************
vector<double> Summary::getDefaults() const {
        vector<double> locations;

        //number of sequences at 2.5%
        const double ptile0_25	= 1+(total * 0.025);
        //number of sequences at 25%
        const double ptile25		= 1+(total * 0.250);
        const double ptile50		= 1+(total * 0.500);
        const double ptile75		= 1+(total * 0.750);
        const double ptile97_5	    = 1+(total * 0.975);
        const double ptile100	    = (total);

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
vector<double> Summary::getValues(map<double, double>& positions) const {
        const vector<double> defaults = getDefaults();
        vector<double> results; results.resize(7,0);
        double meanPosition; meanPosition = 0;
        double totalSoFar = 0;
        double lastValue = 0;

        //minimum
        if ((positions.begin())->first != -1) {
            results[0] = (positions.begin())->first;
        }

        results[1] = results[0]; results[2] = results[0];
        results[3] = results[0]; results[4] = results[0];
        results[5] = results[0];

        for (auto it = positions.begin(); it != positions.end(); it++) {

            float value = it->first; if (isEqual(value, -1)) { value = 0; }
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

        const double meansPosition = meanPosition / static_cast<double>(total);

        // mean
        results.push_back(meansPosition);

        return results;
}
//**************************************************************************
struct seqSumData {
    // "starts" -> (position -> number of seqs with that start position)
    // "ends" -> (position -> number of seqs with that end position)
    map<string, map<double, double> > results;

    vector<float> counts; // abundances for each read, all ones if no count info
    vector<vector<double>> report;
    vector<vector<int>> fastaReport;
    vector<string> columnNames;

    double start;
    double end;

    seqSumData(){
        start = 0;
        end = 0;
    }
    // fasta summary constructor
    seqSumData(const vector<vector<int>>& f, const double st, const double en, const vector<float>& nam,
               const vector<string>& cn) {
        fastaReport = f;
        start = st;
        end = en;
        counts = nam;
        columnNames = cn;
    }

    // contigs and align summary constructor
    seqSumData(const vector<vector<double>>& f, const double st, const double en,
               const vector<float>& nam, const vector<string>& cn) {
        report = f;
        start = st;
        end = en;
        counts = nam;
        columnNames = cn;
    }
};
//**************************************************************************
Rcpp::DataFrame Summary::summarizeFasta(vector<vector<int>> report,
                                  vector<float> counts) {
        numUniques = counts.size();
        total = numUniques;
        if (!counts.empty()) {
            constexpr float initial_sum = 0;
            total = accumulate(counts.begin(), counts.end(), initial_sum);
        }else {
            counts.resize(numUniques, 1);
        }

        createThreadsFasta(report, counts);
        vector<double> nseqs = getDefaults();
        nseqs.push_back(0);

         Rcpp::DataFrame df = Rcpp::DataFrame::create(
             Rcpp::Named("starts") = getValues(results["starts"]),
             Rcpp::_["ends"] = getValues(results["ends"]),
             Rcpp::_["nbases"] = getValues(results["lengths"]),
             Rcpp::_["ambigs"] = getValues(results["ambigs"]),
             Rcpp::_["polymers"] = getValues(results["polymers"]),
             Rcpp::_["numns"] = getValues(results["numns"]),
             Rcpp::_["numseqs"] = nseqs);

         vector<string> rowNames(8, "");
         rowNames[0] = "Minimum:";
         rowNames[1] = "2.5%-tile:";
         rowNames[2] = "25%-tile:";
         rowNames[3] = "Median:   ";
         rowNames[4] = "75%-tile:";
         rowNames[5] = "97.5%-tile:";
         rowNames[6] = "Maximum:";
         rowNames[7] = "Mean:      ";
         df.attr("row.names") = rowNames;
        return (df);
}
//**************************************************************************
// fasta summary data: starts, ends, lengths, ambigs, polymers, numns
void driverSummarize(seqSumData* params) {

    // for each column in dataframe
    for (size_t j = 0; j < params->columnNames.size(); j++) {

        RcppThread::checkUserInterrupt();

        string col = params->columnNames[j];

        // for this processes section of the column
        for (int i = static_cast<int>(params->start); i < static_cast<int>(params->end); i++) {

            const float num = params->counts[i];

            const int thisValue = params->fastaReport[j][i];
            auto it = params->results[col].find(thisValue);
            if (it == params->results[col].end()) {
                params->results[col][thisValue] = num; }
            else { it->second += num; } //add counts
        }
    }
}
//**************************************************************************
void Summary::createThreadsFasta(const vector<vector<int>>& summary,
                                 const vector<float>& counts) {
         //divide reads between processors
         Utils util;
        const vector<pieceOfWork> startEndIndexes = util.divideWork(counts.size(),
                                                              processors);

        vector<RcppThread::Thread*> workerThreads;
        vector<seqSumData*> data;

        vector<string> names(6, "");
        names[0] = "starts"; names[1] = "ends"; names[2] = "lengths";
        names[3] = "ambigs"; names[4] = "polymers"; names[5] = "numns";

        //Launch worker threads
        for (int i = 0; i < processors-1; i++) {

            seqSumData* dataBundle = new seqSumData(summary,
                                                    startEndIndexes[i+1].start,
                                                    startEndIndexes[i+1].end,
                                                    counts, names);
            data.push_back(dataBundle);

            workerThreads.push_back(new RcppThread::Thread(driverSummarize,
                                                           dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(summary, startEndIndexes[0].start,
                                                startEndIndexes[0].end,
                                                counts, names);
        driverSummarize(dataBundle);
        results = dataBundle->results;
        delete dataBundle;

        // for each thread
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();

            // for each column in dataframe
            for (size_t j = 0; j < names.size(); j++) {

                // merge results for this processors section of the column
                for (auto it = data[i]->results[names[j]].begin();
                     it != data[i]->results[names[j]].end(); it++) {
                    auto itMain = results[names[j]].find(it->first);
                    if (itMain == results[names[j]].end()) { //newValue
                        results[names[j]][it->first] = it->second;
                    }else { itMain->second += it->second; } //merge counts
                }
            }

            delete data[i];
            delete workerThreads[i];
        }
}
//**************************************************************************
Rcpp::DataFrame Summary::summarize(Rcpp::DataFrame df, vector<float> counts,
                                   vector<string> dfNames) {

        numUniques = counts.size();
        constexpr float initial_sum = 0;
        total = accumulate(counts.begin(), counts.end(), initial_sum);

        // columns by rows
        vector<vector<double>> report(df.size());
        for (int i = 0; i < df.size(); i++) {
            report[i] = Rcpp::as<std::vector<double> >(df[i]);
        }

        createThreadsReport(report, counts, dfNames);

        Rcpp::DataFrame summaryResult = Rcpp::DataFrame::create();

        for (auto it = results.begin(); it != results.end(); it++){
            summaryResult.push_back(getValues(it->second));
        }
        sort(dfNames.begin(), dfNames.end());
        summaryResult.attr("names") = dfNames;

        vector<string> rowNames(8, "");
        rowNames[0] = "Minimum:";
        rowNames[1] = "2.5%-tile:";
        rowNames[2] = "25%-tile:";
        rowNames[3] = "Median:   ";
        rowNames[4] = "75%-tile:";
        rowNames[5] = "97.5%-tile:";
        rowNames[6] = "Maximum:";
        rowNames[7] = "Mean:      ";
        summaryResult.attr("row.names") = rowNames;

        return (summaryResult);
}
//**************************************************************************
void driverReports(seqSumData* params) {
    // for each column in dataframe
    for (size_t j = 0; j < params->columnNames.size(); j++) {

        RcppThread::checkUserInterrupt();

        string col = params->columnNames[j];

        // for this processes section of the column
        for (size_t i = params->start; i < params->end; i++) {

            const float num = params->counts[i];

            double thisValue = params->report[j][i];
            auto it = params->results[col].find(thisValue);
            if (it == params->results[col].end()) {
                params->results[col][thisValue] = num; }
            else { it->second += num; } //add counts
        }
    }
}
//**************************************************************************
void Summary::createThreadsReport(const vector<vector<double>>& report,
                                   const vector<float>& counts, const vector<string>& names) {
         //divide reads between processors
         Utils util;
        const vector<pieceOfWork> startEndIndexes = util.divideWork(counts.size(),
                                                              processors);

        vector<RcppThread::Thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            seqSumData* dataBundle = new seqSumData(report,
                                                    startEndIndexes[i+1].start,
                                                    startEndIndexes[i+1].end,
                                                    counts, names);
            data.push_back(dataBundle);

            workerThreads.push_back(new RcppThread::Thread(driverReports,
                                                           dataBundle));
        }
        seqSumData* dataBundle = new seqSumData(report,
                                                startEndIndexes[0].start,
                                                startEndIndexes[0].end,
                                                counts, names);

        driverReports(dataBundle);
        results = dataBundle->results;
        delete dataBundle;

        // for each thread
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();

            // for each column in dataframe
            for (size_t j = 0; j < names.size(); j++) {

                // merge results for this processors section of the column
                for (auto it = data[i]->results[names[j]].begin();
                     it != data[i]->results[names[j]].end(); it++) {
                    auto itMain = results[names[j]].find(it->first);
                    if (itMain == results[names[j]].end()) { //newValue
                        results[names[j]][it->first] = it->second;
                    }else { itMain->second += it->second; } //merge counts
                }
            }

            delete data[i];
            delete workerThreads[i];
        }
}
//**************************************************************************
