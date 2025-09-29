#include "phylotree.h"

/******************************************************************************/
PhyloTree::PhyloTree(){
	numNodes = 1;
	numSeqs = 0;
	tree.push_back(TaxNode("Root"));
    tree[0].level = 0;
	maxLevel = 0;

	addSeqToTree("unknown;", 1);
}
/******************************************************************************/
int PhyloTree::addSeqToTree(string seqTaxonomy, float numReps){

    numSeqs += numReps;

    tree[0].total += numReps;
    vector<string> taxons;
    int numLevels = split(seqTaxonomy, ';', back_inserter(taxons));
    util.removeConfidences(taxons);

    if (numLevels > maxLevel) { maxLevel = numLevels; }

    int level = 0;
    int currentNode = 0;

    for(string taxon : taxons) {

        level++;

        auto childPointer = tree[currentNode].children.find(taxon);

        //if the node already exists, move on
        if(childPointer != tree[currentNode].children.end()){
            currentNode = childPointer->second;
            tree[currentNode].total += numReps;
        }else{
            //otherwise, create it
            tree.push_back(TaxNode(taxon));

            tree[currentNode].children[taxon] = numNodes;
            tree[numNodes].level = level;
            tree[numNodes].parent = currentNode;
            tree[numNodes].total = numReps;
            currentNode = numNodes;
            numNodes++;
        }
    }

    return level;
}
/******************************************************************************/
TaxNode PhyloTree::get(int seqIndex){

    if (seqIndex < tree.size()) {
        return tree[seqIndex];
    }
	TaxNode empty;
	return empty;
}
/******************************************************************************/
