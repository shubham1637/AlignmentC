#ifndef AFFINEALIGNMENT_H
#define AFFINEALIGNMENT_H

#include "affinealignobj.h"
#include "alignment.h"
#include <limits>

void doAffineAlignment(float *s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment, AffineAlignObj &affineAlignObj);

AlignedIndices getAffineAlignedIndices(AffineAlignObj &affineAlignObj);

template<class T>
T getOlapAffineAlignStartIndices(T *MatrixM, T *MatrixA, T *MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName);

#endif // AFFINEALIGNMENT_H

