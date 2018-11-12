#include <iostream>
#include <stdio.h>
#include <vector>
#include <limits>

#define NA -1
using namespace std;

void getseqSimMat(string seq1, string seq2, float Match, float MisMatch, float *s){
    int ROW_SIZE = seq1.size();
    int COL_SIZE = seq2.size();
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            seq1[i] == seq2[j] ? *((s+i*COL_SIZE) + j) = Match : *((s+i*COL_SIZE) + j) = MisMatch;
        }
    }
}

void printLengthOfSeq(string seq){
    cout << "The length of sequence is " << seq.size() << endl;
}

template<class T>
void printMatrix(T *s, int ROW_SIZE, int COL_SIZE){
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            cout << *((s+i*COL_SIZE) + j) << " ";
        }
        cout << endl;
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


void getOlapAlignStartIndices1(float *Matrix, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol){
    float maxScore = 0;
    int MaxRowIndex, MaxColIndex;
    for(int i = 0; i < ROW_SIZE; i++){
        if(*((Matrix+i*COL_SIZE) + COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            maxScore = *((Matrix+i*COL_SIZE) + COL_SIZE-1);
        }
    }
    for (int j = 0; j < COL_SIZE; j++){
        if(*((Matrix+(ROW_SIZE-1)*COL_SIZE) + j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            maxScore = *((Matrix+(ROW_SIZE-1)*COL_SIZE) + j);
        }
    }
    OlapStartRow = MaxRowIndex;
    OlapStartCol = MaxColIndex;
}

template<class T>
void getOlapAlignStartIndices(T &Matrix, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol){
    cout << ROW_SIZE << " " << ROW_SIZE << endl;
    float maxScore = 0;
    int MaxRowIndex, MaxColIndex;
    for(int i = 0; i < ROW_SIZE; i++){
        if(Matrix[i][COL_SIZE-1] >= maxScore){
            cout << Matrix[i][COL_SIZE-1]<< endl;
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


const int ROW_SIZE = 3; // If I do not do const, it doesn't pass on to struct and throws error
const int COL_SIZE = 4; // that array bound is not an integer constant before ']'.

struct AlignedIndices{
    vector<int> indexA_aligned;
    vector<int> indexB_aligned;
    vector<float> score;
};

struct AlignObj
{
    vector< vector<char> > Traceback;
    vector< vector<float> > M;
    int signalA_len;
    int signalB_len;
    float Gap;
    bool FreeEndGaps;
};

struct AffineAlignObj
{
    float* M;
    float* A;
    float* B;
    string* Traceback;
    AffineAlignObj(int ROW_SIZE, int COL_SIZE)
    {
        M = new float[ROW_SIZE * COL_SIZE];
        A = new float[ROW_SIZE * COL_SIZE];
        B = new float[ROW_SIZE * COL_SIZE];
        Traceback = new string[3 * ROW_SIZE * COL_SIZE];
    }
    ~AffineAlignObj()
    {
        delete[] M;
        delete[] A;
        delete[] B;
        delete[] Traceback;
    }
    int signalA_len;
    int signalB_len;
    float GapOpen;
    float GapExten;
    bool FreeEndGaps;
};

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

    cout << "M matrix is : " << endl;
    printMatrix((float *)M, signalA_len+1, signalB_len+1);
    printMatrix((char *)Traceback, signalA_len+1, signalB_len+1);

    // Copy cumulative score matrix(M) and traceback matrix into a instance of AlignObj.
    for (int i = 0; i < signalA_len+1; i++) {
        vector<float> rowF; // Create an empty row
        vector<char> rowC; // Create an empty row
        for (int j = 0; j < signalB_len+1; j++) {
            rowF.push_back(M[i][j]); // Add an element (column) to the row
            rowC.push_back(Traceback[i][j]); // Add an element (column) to the row
        }
        alignObj.M.push_back(rowF); // Add the row to the main vector
        alignObj.Traceback.push_back(rowC); // Add the row to the main vector
    }

    return alignObj;
}

AffineAlignObj doAffineAlignment(float *s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment){
    AffineAlignObj affineAlignObj(signalA_len+1, signalB_len+1);
    affineAlignObj.FreeEndGaps = OverlapAlignment;
    affineAlignObj.GapOpen = go;
    affineAlignObj.GapExten = ge;
    affineAlignObj.signalA_len = signalA_len;
    affineAlignObj.signalB_len = signalB_len;

    float M[signalA_len+1][signalB_len+1];
    initializeMatrix(*M, 0, signalA_len+1, signalB_len+1);

    float A[signalA_len+1][signalB_len+1];
    initializeMatrix(*A, 0, signalA_len+1, signalB_len+1);

    float B[signalA_len+1][signalB_len+1];
    initializeMatrix(*B, 0, signalA_len+1, signalB_len+1);

    // enum TbPointer{STOP='S', T='T', D='D', L='L'};
    string Traceback[3][signalA_len+1][signalB_len+1]; // Traceback[0,1,2] = Traceback[TrM,TrA,TrB]
    //initializeMatrix(*Traceback, 'S', signalA_len+1, signalB_len+1);

    // Initialize first row and first column for global and overlap alignment.
    float Inf = std::numeric_limits<float>::infinity();
    for(int i = 0; i<=signalA_len; i++){
        M[i][0] = -Inf;
        B[i][0] = -Inf;
        Traceback[0][i][0] = "SS"; //STOP
        Traceback[2][i][0] = "SS"; //STOP
    }
    for(int j = 0; j<=signalB_len; j++){
        M[0][j] = -Inf;
        A[0][j] = -Inf;
        Traceback[0][0][j] = "SS"; //STOP
        Traceback[1][0][j] = "SS"; //STOP
    }
    M[0][0] = 0;
    if(affineAlignObj.FreeEndGaps == true){
        for(int i = 1; i<=signalA_len; i++){
            A[i][0] = 0;
            Traceback[1][i][0] = "TA"; //TOP A
        }
        for(int j = 1; j<=signalB_len; j++){
            B[0][j] = 0;
            Traceback[2][0][j] = "LB"; //LEFT B
        }
    } else {
        for(int i = 1; i<=signalA_len; i++){
            A[i][0] = -(i-1)*ge - go;
            Traceback[1][i][0] = "TA"; //TOP A
        }
        for(int j = 1; j<=signalB_len; j++){
            B[0][j] = -(j-1)*ge - go;
            Traceback[2][0][j] = "LB"; //LEFT B
        }
    }

    // Perform dynamic programming for affine alignment
    float Diago, gapInA, gapInB;
    for(int j=1; j<=signalB_len; j++ ){
        for(int i=1; i<=signalA_len; i++){
            float sI_1J_1 = *((s+(i-1)*signalB_len) + j-1);
            Diago = M[i-1][j-1] + sI_1J_1;
            gapInA = A[i-1][j-1] + sI_1J_1;
            gapInB = B[i-1][j-1] + sI_1J_1;

            // Calculate recursively for matched alignment
            if(Diago>=gapInA && Diago>=gapInB){
                Traceback[0][i][j] = "DM"; // DM: Diagonal TrM
                M[i][j] = Diago;
            }
            else if (gapInA>=Diago && gapInA>=gapInB){
                Traceback[0][i][j] = "DA"; // DA: Diagonal TrA
                M[i][j] = gapInA;
            }
            else{
                Traceback[0][i][j] = "DB"; // DB: Diagonal TrB
                M[i][j] = gapInB;
            }

            // Calculate recursively for gap in signalB
            if((M[i-1][j]-go) >= (A[i-1][j]-ge) && (M[i-1][j]-go) >= (B[i-1][j]-go)){
                Traceback[1][i][j] = "TM"; // TM: Top TrM
                A[i][j] = M[i-1][j]-go;
            }
            else if ((A[i-1][j]-ge) >= (M[i-1][j]-go) && (A[i-1][j]-ge) >= (B[i-1][j]-go)){
                Traceback[1][i][j] = "TA"; // TA: Top TrA
                A[i][j] = A[i-1][j]-ge;
            }
            else{
                Traceback[1][i][j] = "TB"; // TB: Top TrB
                A[i][j] = B[i-1][j]-go;
            }

            // Calculate recursively for gap in signalA
            if((M[i][j-1]-go) >= (A[i][j-1]-go) && (M[i][j-1]-go) >= (B[i][j-1]-ge)){
                Traceback[2][i][j] = "LM"; // TM: Top TrM
                B[i][j] = M[i][j-1]-go;
            }
            else if ((A[i][j-1]-go) >= (M[i][j-1]-go) && (A[i][j-1]-go) >= (B[i][j-1]-ge)){
                Traceback[2][i][j] = "LA"; // TA: Top TrA
                B[i][j] = A[i][j-1]-go;
            }
            else{
                Traceback[2][i][j] = "LB"; // TB: Top TrB
                B[i][j] = B[i][j-1]-ge;
            }
        }
    }

    cout << "M matrix is : " << endl;
    printMatrix((float *)M, signalA_len+1, signalB_len+1);
    printMatrix((float *)A, signalA_len+1, signalB_len+1);
    printMatrix((float *)B, signalA_len+1, signalB_len+1);
    // printMatrix((string *)Traceback[0], signalA_len+1, signalB_len+1);

    affineAlignObj.M = *M;
    affineAlignObj.A = *A;
    affineAlignObj.B = *B;
    affineAlignObj.Traceback = **Traceback;
    cout << "IN the loop " << endl;
    return affineAlignObj;
}

int main()
{
    float Match=10, MisMatch=-2, go=22, ge=7, gap=go;
    string seq1 = "GCAT";
    string seq2 = "CAGTG";
    int seq1Len = seq1.size();
    int seq2Len = seq2.size();
    float s[seq1Len][seq2Len];
    initializeMatrix(*s, 0, seq1Len, seq2Len);
    getseqSimMat(seq1, seq2, Match, MisMatch, &s[0][0]); // getseqSimMat(seq1, seq2, Match, MisMatch, (float *)s);
    cout << "Similarity matrix is : " << endl;
    printMatrix(*s, seq1Len, seq2Len);
    int signalA_len = sizeof(s) / sizeof(s[0]);
    int signalB_len = sizeof(s[0]) / sizeof(s[0][0]);
    bool OverlapAlignment = true;

    AlignObj alignObj;
    alignObj = doAlignment(*s, seq1Len, seq2Len, gap, OverlapAlignment);
    AffineAlignObj affineAlignObj(seq1Len+1, seq2Len+1); // What if this length is different than used inside the function.
    affineAlignObj = doAffineAlignment(*s, seq1Len, seq2Len, go, ge, OverlapAlignment);
    cout << "Out the loop " << endl;

    cout << "M matrix is : " << endl;
    printMatrix((float *)affineAlignObj.M, signalA_len+1, signalB_len+1);
    printMatrix((float *)affineAlignObj.A, signalA_len+1, signalB_len+1);
    printMatrix((float *)affineAlignObj.B, signalA_len+1, signalB_len+1);

    return 0;
}
