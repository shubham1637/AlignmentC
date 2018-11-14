#include "alignment.h"

AlignedIndices getAlignedIndices(AlignObj alignObj){
    AlignedIndices alignedIdx;
    char TracebackPointer;
    int ROW_IDX = alignObj.signalA_len;
    int COL_IDX = alignObj.signalB_len;

    // Overlap Alignment
    if(alignObj.FreeEndGaps == true){
        getOlapAlignStartIndices(alignObj.M, (alignObj.signalA_len)+1, (alignObj.signalB_len)+1, ROW_IDX, COL_IDX);
        //printMatrix(*M, signalA_len+1, signalB_len+1);
        if(ROW_IDX != alignObj.signalA_len){
            for (int i = alignObj.signalA_len; i>ROW_IDX; i--){
                alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), i);
                alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
                alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[i][COL_IDX]);
            }
        } else if (COL_IDX != alignObj.signalB_len){
            for (int j = alignObj.signalB_len; j>COL_IDX; j--){
                alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
                alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), j);
                alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX][j]);
            }
        }
    }

    TracebackPointer = alignObj.Traceback[ROW_IDX][COL_IDX];
    // Traceback path and align row indices to column indices.
    while(TracebackPointer != 'S'){
        // D: Diagonal, T: Top, L: Left
        if(TracebackPointer == 'D'){
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX][COL_IDX]);
            ROW_IDX = ROW_IDX-1;
            COL_IDX = COL_IDX-1;
        } else if(TracebackPointer == 'T') {
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX][COL_IDX]);
            ROW_IDX = ROW_IDX-1;
        } else {
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), alignObj.M[ROW_IDX][COL_IDX]);
            COL_IDX = COL_IDX-1;
        }
        TracebackPointer = alignObj.Traceback[ROW_IDX][COL_IDX];
        }

    return alignedIdx;
}

AlignObj doAlignment(float *s, int signalA_len, int signalB_len, float gap, bool OverlapAlignment){
    AlignObj alignObj;
    alignObj.FreeEndGaps = OverlapAlignment;
    alignObj.Gap = gap;
    alignObj.signalA_len = signalA_len;
    alignObj.signalB_len = signalB_len;

    float M[signalA_len+1][signalB_len+1];
    initializeMatrix(*M, 0, signalA_len, signalB_len);

    // enum TbPointer{STOP='S', T='T', D='D', L='L'};
    char Traceback[signalA_len+1][signalB_len+1];
    initializeMatrix(*Traceback, 'S', signalA_len+1, signalB_len+1);

    // Initialize first row and first column for global and overlap alignment.
    for(int i = 0; i<=signalA_len; i++){
        M[i][0] = -i*gap;
        Traceback[i][0] = 'T'; //Top
    }
    for(int j = 0; j<=signalB_len; j++){
        M[0][j] = -j*gap;
        Traceback[0][j] = 'L'; //Left
    }
    Traceback[0][0] = 'S'; //STOP

    // Perform dynamic programming for alignment
    float Diago, gapInA, gapInB;
    for(int i=1; i<=signalA_len; i++ ){
        for(int j=1; j<=signalB_len; j++){
            Diago = M[i-1][j-1] + *((s+(i-1)*signalB_len) + j-1);
            gapInA = M[i-1][j] - gap;
            gapInB = M[i][j-1] - gap;
            if(Diago>=gapInA && Diago>=gapInB){
                Traceback[i][j] = 'D'; // D: Diagonal
                M[i][j] = Diago;
            }
            else if (gapInA>=Diago && gapInA>=gapInB){
                Traceback[i][j] = 'L'; // L: Left
                M[i][j] = gapInA;
            }
            else{
                Traceback[i][j] = 'T'; // T: Top
                M[i][j] == gapInB;
            }
        }
    }

    std::cout << "M matrix is : " << std::endl;
    printMatrix((float *)M, signalA_len+1, signalB_len+1);
    printMatrix((char *)Traceback, signalA_len+1, signalB_len+1);

    // Copy cumulative score matrix(M) and traceback matrix into a instance of AlignObj.
    for (int i = 0; i < signalA_len+1; i++) {
        std::vector<float> rowF; // Create an empty row
        std::vector<char> rowC; // Create an empty row
        for (int j = 0; j < signalB_len+1; j++) {
            rowF.push_back(M[i][j]); // Add an element (column) to the row
            rowC.push_back(Traceback[i][j]); // Add an element (column) to the row
        }
        alignObj.M.push_back(rowF); // Add the row to the main vector
        alignObj.Traceback.push_back(rowC); // Add the row to the main vector
    }

    return alignObj;
}

template<class T>
void printMatrix(T *s, int ROW_SIZE, int COL_SIZE){
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            std::cout << *((s+i*COL_SIZE) + j) << " ";
        }
        std::cout << std::endl;
    }
}

template<class T1, class T2> // Don't know why do I have to specify this line everytime. It is just the rule.
void initializeMatrix(T1 *s, T2 val, int ROW_SIZE, int COL_SIZE){
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            *((s+i*COL_SIZE) + j) = val;
        }
    }
}

template<class T>
void getOlapAlignStartIndices(T &Matrix, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol){
    std::cout << ROW_SIZE << " " << ROW_SIZE << std::endl;
    float maxScore = 0;
    int MaxRowIndex, MaxColIndex;
    for(int i = 0; i < ROW_SIZE; i++){
        if(Matrix[i][COL_SIZE-1] >= maxScore){
            std::cout << Matrix[i][COL_SIZE-1]<< std::endl;
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            maxScore = Matrix[i][COL_SIZE-1];
        }
    }
    for (int j = 0; j < COL_SIZE; j++){
        if(Matrix[ROW_SIZE-1][j] >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            maxScore = Matrix[ROW_SIZE-1][j];
        }
    }
    OlapStartRow = MaxRowIndex;
    OlapStartCol = MaxColIndex;
}
