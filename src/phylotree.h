#ifndef DOTAXONOMY_H
#define DOTAXONOMY_H

#include "utils.h"

/**************************************************************************************************/

struct TaxNode {
	//vector<string> accessions;	//names of seqs in this branch of tree
	map<string, int> children;  //childs name to index in tree

	int parent, childNumber, level, total;
	string name;

	TaxNode(string n) : name(n) {
	    initialize();
	}
	TaxNode() : name("") { initialize(); }

	void print(ostream& out) {
	    out << "name / total / level " << name << '\t' << total << '\t' << level << endl;
	    out << "parent " << parent << endl;
	    out << "children " << toString(getKeys(children), ',') << endl;
	}

	void initialize() {
	    level = 0;
	    parent = -1;
        total = 0;
	}
};

/**************************************************************************************************/

class PhyloTree {


public:
	PhyloTree();
	~PhyloTree() = default;

	int addSeqToTree(string, float);
	TaxNode get(int seqIndex);

	TaxNode getRoot()       {   return tree[0];     }
	int getMaxLevel()		{	return maxLevel;	}
	float getNumSeqs()		{	return numSeqs;		}
	int getNumNodes()		{	return (int)tree.size();	}

private:

    Utils util;
	vector<TaxNode> tree;

 	int numNodes;
	float numSeqs;
	int maxLevel;
};

/**************************************************************************************************/

#endif


