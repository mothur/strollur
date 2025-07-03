#ifndef utils_h
#define utils_h

#include "dataset.h"

class Utils {

public:

    Utils() = default;
    ~Utils() = default;

    vector<pieceOfWork> divideWork(double numItems, int& numProcessors) {
        // divide work between processors
        vector<pieceOfWork> work;

        if (numItems < numProcessors) { numProcessors = numItems; }
        size_t startIndex = 0;

        for (size_t remainingProcessors = numProcessors; remainingProcessors > 0;
        remainingProcessors--) {

            //case for last processor
            size_t numToProcess = numItems;
            if (remainingProcessors != 1) {
                numToProcess = ceil(numItems / remainingProcessors);
            }
            work.push_back(pieceOfWork(startIndex, (startIndex+numToProcess)));
            startIndex += numToProcess;
            numItems -= numToProcess;
        }
        return work;
    }

    bool isPositiveNumeric(string s){
        bool numeric = false;

        if (s == "") { numeric = false;  }
        else if(s.find_first_not_of("0123456789.") == string::npos) {
            numeric = true;
        }

        return numeric;
    }

    float round2Places(float var) {
        // 37.66666 * 100 =3766.66
        // 3766.66 + .5 =3767.16    for rounding off value
        // then type cast to int so value is 3767
        // then divided by 100 so the value converted into 37.67
        float value = (int)(var * 100 + .5);
        return (float)value / 100;
    }

    // taxonomy utils
    int removeConfidence(string& taxon) {
        int confidence = 0;
        int pos = taxon.find_last_of('(');
        if (pos != -1) {
            //is it a number?
            int pos2 = taxon.find_last_of(')');
            if (pos2 != -1) {
                string temp = taxon.substr(pos+1, (pos2-(pos+1)));
                if (isPositiveNumeric(temp)) {
                    taxon = taxon.substr(0, pos); //rip off confidence
                    convert(temp, confidence);
                }
            }
        }
        return confidence;
    }
    vector<int> removeConfidences(vector<string>& taxons) {
        vector<int> confidences(taxons.size(), 0);

        for (int i = 0; i < taxons.size(); i++) {
            confidences[i] = removeConfidence(taxons[i]);
        }
        return confidences;
    }

    void addUnclassifieds(vector<string>& taxons,
                          vector<int>& confidences, int maxlevel) {

        if (taxons.size() == maxlevel) { return; }

        int index = 0;
        // find last taxon, end early for unclassifieds
        for (int i = 0; i < taxons.size(); i++) {
            index = i;
            if (taxons[i] == "unclassified"){ index--; break; }
        }
        int level = index+1;
        string cTax = taxons[index] + "_unclassified";

        //add "unclassified" until you reach maxLevel
        while (level < maxlevel) {
            taxons.push_back(cTax);
            level++;
        }

        if (!confidences.empty()) {
            int level = index+1;
            int cConfidence = confidences[index];

            while (level < maxlevel) {
                confidences.push_back(cConfidence);
                level++;
            }
        }
    }

private:


};

#endif /* utils_h */

