#ifndef AFFINEALIGNMENT_H
#define AFFINEALIGNMENT_H

#include "affinealignobj.h"
#include "alignment.h"
#include <limits>

// TODO: needs more comments!
// TODO: the similarity matrix should also be an object that knows about its size (n_rows, n_columns)
// TODO: why not return an affine alignment object? instead of making the fxn void
void doAffineAlignment(float *s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment, AffineAlignObj &affineAlignObj);

AlignedIndices getAffineAlignedIndices(AffineAlignObj &affineAlignObj);

// TODO: is this an internal function? if so, move to cpp
template<class T>
T getOlapAffineAlignStartIndices(T *MatrixM, T *MatrixA, T *MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName);

#endif // AFFINEALIGNMENT_H

