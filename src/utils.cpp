
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
bool Utils::isPositiveNumeric(const string& s){
    bool numeric = false;

    if (s.empty()) { numeric = false;  }
    else if(s.find_first_not_of("0123456789.") == string::npos) {
        numeric = true;
    }

    return numeric;
}
/******************************************************************************/
// taxonomy utils
int Utils::removeConfidence(string& taxon) {
    int confidence = 0;
    const int pos = taxon.find_last_of('(');
    if (pos != -1) {
        //is it a number?
        const int pos2 = taxon.find_last_of(')');
        if (pos2 != -1) {
            const string temp = taxon.substr(pos+1, (pos2-(pos+1)));
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

    for (size_t i = 0; i < taxons.size(); i++) {
        confidences[i] = removeConfidence(taxons[i]);
    }
    return confidences;
}
// /******************************************************************************/
// void Utils::addUnclassifieds(string& taxon, int maxlevel) {
//     vector<string> taxons;
//     split(taxon, ';', back_inserter(taxons));
//     vector<int> confidences = removeConfidences(taxons);
//     addUnclassifieds(taxons, confidences, maxlevel);
//
//     taxon = "";
//     for (int i = 0; i < maxlevel; i++) {
//         taxon += taxons[i] + "(" + toString(confidences[i]) + ");";
//     }
// }
/******************************************************************************/
void Utils::addUnclassifieds(vector<string>& taxons,
                      vector<int>& confidences, const int maxlevel) {

    if (static_cast<int>(taxons.size()) == maxlevel) { return; }

    int index = 0;
    // find last taxon, end early for unclassifieds
    for (int i = 0; i < static_cast<int>(taxons.size()); i++) {
        index = i;
        if (taxons[i] == "unclassified"){ index--; break; }
    }
    int level = index+1;
    const string cTax = taxons[index] + "_unclassified";

    //add "unclassified" until you reach maxLevel
    while (level < maxlevel) {
        taxons.push_back(cTax);
        level++;
    }

    if (!confidences.empty()) {
        int level = index+1;
        const int cConfidence = confidences[index];

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
            for (size_t i = 0; i < stax.size(); i++) {
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

    for (size_t i = 0; i < tax.size(); i++) {
        tax[i] = removeQuotes(tax[i]);
    }
}
/**************************************************************************************************/
string Utils::removeQuotes(const string &tax) {
    string taxon;
    string newTax = "";

    for (size_t i = 0; i < tax.length(); i++) {
        if ((tax[i] != '\'') && (tax[i] != '\"')) { newTax += tax[i]; }
    }

    return newTax;
}
/**************************************************************************************************/
bool Utils::searchTax(const vector<string> &userTaxons,
                      const vector<int> &userConfidences,
                      const vector<bool>& taxonsHasConfidence,
                      const vector< vector<string> > &searchTaxons,
                      const vector< vector<int> > &searchConfidenceThresholds) {

    const bool userDataHasConfidence = (sum(userConfidences) != 0);

    for (size_t j = 0; j < searchTaxons.size(); j++) {

        const bool foundTaxonMatch = findTaxon(userTaxons, searchTaxons[j]);
        const int userTaxonSize = static_cast<int>(userTaxons.size());
        if (foundTaxonMatch) {
            // searchTaxon or user taxons don't include confidence scores so
            // ignore them
            if (!taxonsHasConfidence[j] || !userDataHasConfidence) {
                // since you belong to at least one of the taxons we want you
                // are included so no need to search for other
                return true;
            }else {
                bool found = false;

                // the usersTaxon is most likely longer than the searchTaxons,
                // and searchTaxon[0] may relate to userTaxon[4]
                // we want to "line them up", so we will find the the index
                // where the searchstring starts
                int index = 0;
                for (int i = 0; i < userTaxonSize; i++) {

                    if (userTaxons[i] == searchTaxons[j][0]) {
                        index = i;
                        int spot = 0;
                        bool goodspot = true;

                        // is this really the start, or are we dealing with a
                        // taxon of the same name?
                        while ((spot < static_cast<int>(searchTaxons[j].size())) &&
                               ((i+spot) < userTaxonSize)) {
                            if (userTaxons[i+spot] != searchTaxons[j][spot]) {
                                goodspot = false; break;
                            }else { spot++; }
                        }

                        if (goodspot) { break; }
                    }
                }

                for (int i = 0; i < static_cast<int>(searchTaxons[j].size()); i++) {

                    if ((i+index) < userTaxonSize) {
                        //is the users cutoff less than the search taxons
                        if (userConfidences[i+index] < searchConfidenceThresholds[j][i]) {
                            found = true;
                            break;
                        }
                    }else { found = true; break; }
                }

                //passed the test so add you
                if (found) { return true; }
            }
        }
    }

    return false;
}
/**************************************************************************************************/
