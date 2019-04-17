#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <stdio.h>
#include <vector>

#include "affinealignobj.h"
#define NA 0

struct AlignedIndices{
    std::vector<int> indexA_aligned;
    std::vector<int> indexB_aligned;
    std::vector<float> score;
};

struct AlignObj
{
    std::vector< std::vector<char> > Traceback;
    std::vector< std::vector<float> > M;
    int signalA_len;
    int signalB_len;
    float Gap;
    bool FreeEndGaps;
};

AlignedIndices getAlignedIndices(AlignObj alignObj);

/** @brief Perform (simple) alignment on similarly matrix s
 *
 * @param s Input similarly matrix
 * @param signalA_len The number of rows
 * @param signalB_len The number of columns
 * @param gap The gap penalty
 * @param OverlapAlignment Whether end gaps should be free
 *
 * @return AlignObj with filled out matrix M and Traceback
 */
AlignObj doAlignment(float *s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment);

template<class T>
void getOlapAlignStartIndices(T &Matrix, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol);

/** @brief Initializes matrix with a given value
 *
 * @param T1 A pointer to the matrix to initialize
 * @param T2 The value to use for initialization
 * @param ROW_SIZE The number of rows
 * @param COL_SIZE The number of columns
 */
template<class T1, class T2> // Don't know why do I have to specify this line everytime. It is just the rule.
void initializeMatrix(T1 *s, T2 val, int ROW_SIZE, int COL_SIZE);

template<class T>
void printMatrix(T *s, int ROW_SIZE, int COL_SIZE){
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            std::cout << *((s+i*COL_SIZE) + j) << " ";
        }
        std::cout << std::endl;
    }
}

#endif // ALIGNMENT_H
