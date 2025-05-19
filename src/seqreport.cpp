

#include "seqreport.h"

/******************************************************************************/
vector<int> SeqReport::getReport(string seq) {
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
/******************************************************************************/
vector<vector<int>> SeqReport::getReport(vector<string>& seqs) {
    vector<vector<int>> results(6, vector<int>(seqs.size()));

    for (int i = 0; i < seqs.size(); i++) {
        vector<int> thisSeqsResults = getReport(seqs[i]);
        results[0][i] = thisSeqsResults[0];
        results[1][i] = thisSeqsResults[1];
        results[2][i] = thisSeqsResults[2];
        results[3][i] = thisSeqsResults[3];
        results[4][i] = thisSeqsResults[4];
        results[5][i] = thisSeqsResults[5];
    }

    return results;
}
/******************************************************************************/
// assumes starts, ends, numbases, ambigs, polymers and numns
// are all the same size
void SeqReport::addReports(vector<string>& seqs, vector<int>& starts,
                           vector<int>& ends, vector<int>& numbases,
                           vector<int>& ambigs, vector<int>& polymers,
                           vector<int>& numns){

    int originalSize = starts.size();
    int numSeqs = seqs.size();

    starts.resize(originalSize+numSeqs);
    ends.resize(originalSize+numSeqs);
    numbases.resize(originalSize+numSeqs);
    ambigs.resize(originalSize+numSeqs);
    polymers.resize(originalSize+numSeqs);
    numns.resize(originalSize+numSeqs);

    for (int i = originalSize; i < originalSize+numSeqs; i++) {
        vector<int> thisSeqsResults = getReport(seqs[i]);
        starts[i] = thisSeqsResults[0];
        ends[i] = thisSeqsResults[1];
        numbases[i] = thisSeqsResults[2];
        ambigs[i] = thisSeqsResults[3];
        polymers[i] = thisSeqsResults[4];
        numns[i] = thisSeqsResults[5];
    }
}
/******************************************************************************/
int SeqReport::getStart(string seq) {
    int startPos = 1;
    for(int i = 0; i < seq.length(); i++) {
        if(!isgap(seq[i])){
            startPos = i + 1;
            break;
        }
    }
    return startPos;
}
/******************************************************************************/
int SeqReport::getEnd(string seq) {
    int endPos = 1;
    for(int i = seq.length()-1; i >= 0; i--){
        if(!isgap(seq[i])){
            endPos = i + 1;
            break;
        }
    }
    return endPos;
}
/******************************************************************************/
int SeqReport::getNumAmbigs(string seq) {
    int ambigBases = 0;
    for(int i = 0; i < seq.length(); i++) {
        if(!isgap(seq[i]) && !isacgt(seq[i])){
            ambigBases++;
        }
    }
    return ambigBases;
}
/******************************************************************************/
int SeqReport::getLongestHomopolymer(string seq) {
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
/******************************************************************************/
int SeqReport::getNumbases(string seq) {
    int numBases = 0;
    for(int i = 0; i < seq.length(); i++) {
        if(!isgap(seq[i])){
            numBases++;
        }
    }
    return numBases;
}
/******************************************************************************/
int SeqReport::getNumns(string seq) {
    int numNs = 0;
    for (int i = 0; i < seq.length(); i++) {
        if(seq[i] == 'N') { numNs++; }
    }
    return numNs;
}
/******************************************************************************/
