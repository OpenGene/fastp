#include "nucleotidetree.h"

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
    //cout << base;
    //cout << count;
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
        char base = seq[i] & 0x07;
        if(curNode->children[base] == NULL) {
            curNode->children[base] = new NucleotideNode();
            curNode->children[base]->base = seq[i];
        }
        curNode->children[base]->count++;
        curNode = curNode->children[base];
    }
}

bool NucleotideTree::test() {
    NucleotideTree tree(NULL);
    tree.addSeq("AAAATTTTCCCCGGGG");
    tree.addSeq("AAAATTTTGGGGCCCC");
    tree.addSeq("AAAATTTTGGGGCCCT");
    tree.mRoot->dfs();
    return true;
}