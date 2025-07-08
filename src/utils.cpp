
#include "utils.h"

/******************************************************************************/
vector<pieceOfWork> Utils::divideWork(double numItems, int& numProcessors) {
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
/******************************************************************************/
bool Utils::isPositiveNumeric(string s){
    bool numeric = false;

    if (s == "") { numeric = false;  }
    else if(s.find_first_not_of("0123456789.") == string::npos) {
        numeric = true;
    }

    return numeric;
}
/******************************************************************************/
float Utils::round2Places(float var) {
    // 37.66666 * 100 =3766.66
    // 3766.66 + .5 =3767.16    for rounding off value
    // then type cast to int so value is 3767
    // then divided by 100 so the value converted into 37.67
    float value = (int)(var * 100 + .5);
    return (float)value / 100;
}
/******************************************************************************/
// taxonomy utils
int Utils::removeConfidence(string& taxon) {
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
/******************************************************************************/
vector<int> Utils::removeConfidences(vector<string>& taxons) {
    vector<int> confidences(taxons.size(), 0);

    for (int i = 0; i < taxons.size(); i++) {
        confidences[i] = removeConfidence(taxons[i]);
    }
    return confidences;
}
/******************************************************************************/
void Utils::addUnclassifieds(string& taxon, int maxlevel) {
    vector<string> taxons;
    split(taxon, ';', back_inserter(taxons));
    vector<int> confidences = removeConfidences(taxons);
    addUnclassifieds(taxons, confidences, maxlevel);

    taxon = "";
    for (int i = 0; i < maxlevel; i++) {
        taxon += taxons[i] + "(" + toString(confidences[i]) + ");";
    }
}
/******************************************************************************/
void Utils::addUnclassifieds(vector<string>& taxons,
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

/******************************************************************************/
/*M02352_41_000000000-AT06G_1_2104_18738_21630 Eukaryota(100);Archaeplastida(100);Chloroplastida(100);Chlorophyta(100);Mamiellophyceae(100);Mamiellales(100);Ostreococcus(100);Ostreococcus tauri(100);

 When I run remove.lineage with:
 taxon=Chloroplast-Mitochondria-unknown-Bacteria-Archaea-Metazoa-Charophyta

 The word "Chloroplast" in the taxon string gets matched to the lineage Chloroplastida in the taxonomy (above) and wipes out all of the green algae.*/

bool Utils::findTaxon(vector<string> tax, vector<string> stax) {

    removeQuotes(tax); removeQuotes(stax);

    //looking to find something like "unknown" or "Proteobacteria"
    if (stax.size() == 1) {
        string searchTax = stax[0];
        auto it = find_if(tax.begin(), tax.end(),
                          [&searchTax](const string& obj) {
                              return obj == searchTax;
                          });

        if (it != tax.end()) { return true; }
        else { return false; }

    }else { //looking to find something like
        // "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Anaplasmataceae;Wolbachia;"

        //we are looking for a more specific taxonomy, not a match
        if (stax.size() > tax.size()) { return false; } //we are looking for a more specific taxonomy, not a match
        else {
            for (int i = 0; i < stax.size(); i++) {
                if (stax[i] != tax[i]) { return false; }
            }
            return true;
        }
    }

    return false;
}
/**************************************************************************************************/
void Utils::removeQuotes(vector<string>& tax) {
    string taxon;
    string newTax = "";

    for (int i = 0; i < tax.size(); i++) {
        tax[i] = removeQuotes(tax[i]);
    }
}
/**************************************************************************************************/
string Utils::removeQuotes(string tax) {
    string taxon;
    string newTax = "";

    for (int i = 0; i < tax.length(); i++) {
        if ((tax[i] != '\'') && (tax[i] != '\"')) { newTax += tax[i]; }
    }

    return newTax;
}
/**************************************************************************************************/
bool Utils::searchTax(vector<string> userTaxons,
                      vector<int> userConfidences,
                      vector<bool> taxonsHasConfidence,
                      vector< vector<string> > searchTaxons,
                      vector< vector<int> > searchConfidenceThresholds) {

    bool userDataHasConfidence = (sum(userConfidences) != 0);

    for (int j = 0; j < searchTaxons.size(); j++) {

        bool foundTaxonMatch = findTaxon(userTaxons, searchTaxons[j]);

        if (foundTaxonMatch) {
            // searchTaxon or user taxons don't include confidence scores so
            // ignore them
            if (!taxonsHasConfidence[j] || !userDataHasConfidence) {
                // since you belong to at least one of the taxons we want you
                // are included so no need to search for other
                return true;
            }else {
                bool good = true;

                // the usersTaxon is most likely longer than the searchTaxons,
                // and searchTaxon[0] may relate to userTaxon[4]
                // we want to "line them up", so we will find the the index
                // where the searchstring starts
                int index = 0;
                for (int i = 0; i < userTaxons.size(); i++) {

                    if (userTaxons[i] == searchTaxons[j][0]) {
                        index = i;
                        int spot = 0;
                        bool goodspot = true;

                        // is this really the start, or are we dealing with a
                        // taxon of the same name?
                        while ((spot < searchTaxons[j].size()) &&
                               ((i+spot) < userTaxons.size())) {
                            if (userTaxons[i+spot] != searchTaxons[j][spot]) {
                                goodspot = false; break;
                            }else { spot++; }
                        }

                        if (goodspot) { break; }
                    }
                }

                for (int i = 0; i < searchTaxons[j].size(); i++) {

                    if ((i+index) < userTaxons.size()) {
                        //is the users cutoff less than the search taxons
                        if (userConfidences[i+index] < searchConfidenceThresholds[j][i]) {
                            good = false;
                            break;
                        }
                    }else { good = false; break; }
                }

                //passed the test so add you
                if (good) { return true; }
            }
        }
    }

    return false;
}
/**************************************************************************************************/
