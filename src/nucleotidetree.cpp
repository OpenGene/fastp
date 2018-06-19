#include "nucleotidetree.h"
#include <sstream>

NucleotideNode::NucleotideNode(){
    count = 0;
    base = 'N';
    memset(children, 0, sizeof(NucleotideNode*)*8);
}
NucleotideNode::~NucleotideNode(){
    for(int i=0; i<8; i++) {
        if(children[i])
            delete children[i];
    }
}
void NucleotideNode::dfs() {
    //cerr << base;
    //cerr << count;
    printf("%c", base);
    printf("%d", count);
    bool hasChild = false;
    for(int i=0; i<8; i++) {
        if(children[i]){
            children[i]->dfs();
            hasChild = true;
        }
    }
    if(!hasChild) {
        printf("\n");
    }
}

NucleotideTree::NucleotideTree(Options* opt){
    mOptions = opt;
    mRoot = new NucleotideNode();
}


NucleotideTree::~NucleotideTree(){
    delete mRoot;
}

void NucleotideTree::addSeq(string seq) {
    NucleotideNode* curNode = mRoot;
    for(int i=0; i<seq.length(); i++) {
        if(seq[i] == 'N')
            break;
        char base = seq[i] & 0x07;
        if(curNode->children[base] == NULL) {
            curNode->children[base] = new NucleotideNode();
            curNode->children[base]->base = seq[i];
        }
        curNode->children[base]->count++;
        curNode = curNode->children[base];
    }
}

string NucleotideTree::getDominantPath(bool& reachedLeaf) {
    stringstream ss;
    const double RATIO_THRESHOLD = 0.95;
    const int NUM_THRESHOLD = 50;
    NucleotideNode* curNode = mRoot;
    while(true) {
        int total = 0;
        for(int i=0; i<8; i++) {
            if(curNode->children[i] != NULL)
                total += curNode->children[i]->count;
        }
        if(total < NUM_THRESHOLD)
            break;
        bool hasDominant = false;
        for(int i=0; i<8; i++) {
            if(curNode->children[i] == NULL)
                continue;
            if(curNode->children[i]->count / (double)total >= RATIO_THRESHOLD) {
                hasDominant = true;
                ss << curNode->children[i]->base;
                curNode = curNode->children[i];
                break;
            }
        }
        if(!hasDominant) {
            reachedLeaf = false;
            break;
        }
    }
    return ss.str();

}

bool NucleotideTree::test() {
    NucleotideTree tree(NULL);
    for(int i=0; i<100; i++) {
        tree.addSeq("AAAATTTT");
        tree.addSeq("AAAATTTTGGGG");
        tree.addSeq("AAAATTTTGGGGCCCC");
        tree.addSeq("AAAATTTTGGGGCCAA");
    }
    tree.addSeq("AAAATTTTGGGACCCC");

    bool reachedLeaf = true;
    string path = tree.getDominantPath(reachedLeaf);
    printf("%s\n", path.c_str());
    return path == "AAAATTTTGGGGCC";
}