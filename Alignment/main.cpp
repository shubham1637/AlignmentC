#include <iostream>
#include <stdio.h>
#include <vector>
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


int main()
{
    int m[ROW_SIZE][COL_SIZE]=
        {
            {1,2,3,4},
            {5,6,7,8},
            {13,14,15,16}
        };
    int VecMat[ROW_SIZE*COL_SIZE]= {1,5,13, 2,6,14, 3,7,15, 4,8,16};

    int k[ROW_SIZE][COL_SIZE] = {0};
    int ROW_IDX;
    int COL_IDX;
    for(int i=0; i < ROW_SIZE*COL_SIZE; i++){
        ROW_IDX = i%ROW_SIZE;
        COL_IDX = i/ROW_SIZE;
        k[ROW_IDX][COL_IDX] = 2*VecMat[i];
    }
    cout << "Printing K" << endl;
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            cout << k[i][j] << " ";
        }
        cout << endl;
    }

    //NeedleAlignObj doGlobalAlignment(float s, float gap){
    //}
    float Match=10, MisMatch=-2, go=22, ge=7, gap=go;
    string seq1 = "GCAT";
    string seq2 = "CAGTG";
    int seq1Len = seq1.size();
    int seq2Len = seq2.size();
    float s[seq1Len][seq2Len];
    initializeMatrix(*s, 0, seq1Len, seq2Len);
    printMatrix(*s, seq1Len, seq2Len);
    getseqSimMat(seq1, seq2, Match, MisMatch, &s[0][0]); // getseqSimMat(seq1, seq2, Match, MisMatch, (float *)s);
    cout << "Similarity matrix is : " << endl;
    printMatrix(*s, seq1Len, seq2Len);
    cout << sizeof(s) / sizeof(s[0]) << endl;
    cout << sizeof(s[0]) / sizeof(s[0][0]) << endl;
    bool OverlapAlignment = true;

    AlignObj alignObj;
    alignObj = doAlignment(*s, seq1Len, seq2Len, gap, OverlapAlignment);

    AlignedIndices alignedIdx;
    alignedIdx = getAlignedIndices(alignObj);

    int ChromA_Len = sizeof(s) / sizeof(s[0]);
    int ChromB_Len = sizeof(s[0]) / sizeof(s[0][0]);

    cout << endl;
    for (vector<float>::iterator it = alignedIdx.score.begin(); it != alignedIdx.score.end(); it++)
        cout << *it << " ";

    cout << endl;
    for (vector<int>::iterator it = alignedIdx.indexA_aligned.begin(); it != alignedIdx.indexA_aligned.end(); it++)
        cout << *it << " ";

    cout << endl;
    for (vector<int>::iterator it = alignedIdx.indexB_aligned.begin(); it != alignedIdx.indexB_aligned.end(); it++)
        cout << *it << " ";
    cout << endl;

    return 0;
}
