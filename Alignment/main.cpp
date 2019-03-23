#include <iostream>
#include <stdio.h>
#include <vector>

#include "alignment.h"
#include "affinealignobj.h"
#include "affinealignment.h"

void getseqSimMat(std::string seq1, std::string seq2, float Match, float MisMatch, float *s){
    int ROW_SIZE = seq1.size();
    int COL_SIZE = seq2.size();
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            seq1[i] == seq2[j] ? *((s+i*COL_SIZE) + j) = Match : *((s+i*COL_SIZE) + j) = MisMatch;
        }
    }
}

void printLengthOfSeq(std::string seq){
    std::cout << "The length of sequence is " << seq.size() << std::endl;
}

int main()
{
    float Match=10, MisMatch=-2, go=22, ge=7, gap=go;
    std::string seq1 = "GCAT";
    std::string seq2 = "CAGTG";
    int seq1Len = seq1.size();
    int seq2Len = seq2.size();
    float s[seq1Len][seq2Len];
    initializeMatrix(*s, 0, seq1Len, seq2Len);
    getseqSimMat(seq1, seq2, Match, MisMatch, &s[0][0]); // getseqSimMat(seq1, seq2, Match, MisMatch, (float *)s);
    std::cout << "Similarity matrix is : " << std::endl;
    printMatrix(*s, seq1Len, seq2Len);
    int signalA_len = sizeof(s) / sizeof(s[0]);
    int signalB_len = sizeof(s[0]) / sizeof(s[0][0]);
    bool OverlapAlignment = true;

    AlignObj alignObj;
    alignObj = doAlignment(*s, seq1Len, seq2Len, gap, OverlapAlignment);
    AffineAlignObj affineAlignObj(seq1Len+1, seq2Len+1); // What if this length is different than used inside the function.
    AffineAlignObjNew affineAlignObjNew(seq1Len+1, seq2Len+1);
    doAffineAlignment(*s, seq1Len, seq2Len, go, ge, OverlapAlignment, affineAlignObj);
    std::cout << "Out the loop " << std::endl;

//    AffineAlignObj affineAlignObj2(seq1Len+1, seq2Len+1);
//    affineAlignObj2 = affineAlignObj;
//    doAffineAlignment(*s, seq1Len, seq2Len, 3, 3, OverlapAlignment, affineAlignObj);

    std::cout << "M matrix is : " << std::endl;
    printMatrix((float *)affineAlignObj.M, signalA_len+1, signalB_len+1);
    printMatrix((float *)affineAlignObj.A, signalA_len+1, signalB_len+1);
    printMatrix((float *)affineAlignObj.B, signalA_len+1, signalB_len+1);
    printMatrix((TracebackType *)affineAlignObj.Traceback, signalA_len+1, signalB_len+1);
    printMatrix((TracebackType *)affineAlignObj.Traceback+1*((signalA_len+1)*(signalB_len+1)), signalA_len+1, signalB_len+1);
    printMatrix((TracebackType *)affineAlignObj.Traceback+2*((signalA_len+1)*(signalB_len+1)), signalA_len+1, signalB_len+1);

//    std::cout << "M matrix is : " << std::endl;
//    printMatrix((float *)affineAlignObj2.M, signalA_len+1, signalB_len+1);
//    printMatrix((float *)affineAlignObj2.A, signalA_len+1, signalB_len+1);
//    printMatrix((float *)affineAlignObj2.B, signalA_len+1, signalB_len+1);

    AlignedIndices alignedIdx;
    alignedIdx = getAffineAlignedIndices(affineAlignObj);

    std::cout << std::endl;
    for (std::vector<float>::iterator it = alignedIdx.score.begin(); it != alignedIdx.score.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;
    for (std::vector<int>::iterator it = alignedIdx.indexA_aligned.begin(); it != alignedIdx.indexA_aligned.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;
    for (std::vector<int>::iterator it = alignedIdx.indexB_aligned.begin(); it != alignedIdx.indexB_aligned.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;

    return 0;
}
