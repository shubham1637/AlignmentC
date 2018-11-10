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


void getOlapAlignStartIndices(float *Matrix, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol){
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

int main()
{
    const int ROW_SIZE = 3; // If I do not do const, it doesn't pass on to struct and throws error
    const int COL_SIZE = 4; // that array bound is not an integer constant before ']'.

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

    struct AlignedIndices{
        vector<int> indexA_aligned;
        vector<int> indexB_aligned;
        vector<float> score;
    };

    struct NeedleAlignObj
    {
        AlignedIndices Traceback[1];
        int M[ROW_SIZE][COL_SIZE];
        float Gap;
    };

    struct AlignObj
    {
        AlignedIndices Traceback[1];
        float M[ROW_SIZE][COL_SIZE];
        float Gap;
        bool FreeEndGaps;
    };


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


    int ChromA_Len = sizeof(s) / sizeof(s[0]);
    int ChromB_Len = sizeof(s[0]) / sizeof(s[0][0]);
    float M[ChromA_Len+1][ChromB_Len+1];
    initializeMatrix(*M, 0, seq1Len, seq2Len);
    cout << "M matrix is : " << endl;
    printMatrix((float *)M, ChromA_Len+1, ChromB_Len+1);

    // enum TbPointer{STOP='S', T='T', D='D', L='L'};
    char Traceback[ChromA_Len+1][ChromB_Len+1];
    initializeMatrix(*Traceback, 'S', ChromA_Len+1, ChromB_Len+1);
    cout <<"Traceback matrix" <<endl;
    printMatrix((char *)Traceback, ChromA_Len+1, ChromB_Len+1);

    for(int i = 0; i<=ChromA_Len; i++){
        M[i][0] = -i*gap;
        Traceback[i][0] = 'T';
    }
    for(int j = 0; j<=ChromB_Len; j++){
        M[0][j] = -j*gap;
        Traceback[0][j] = 'L';
    }
    Traceback[0][0] = 'S';
    //printMatrix(*M, ChromA_Len+1, ChromB_Len+1);
    //printMatrix((char *)Traceback, ChromA_Len+1, ChromB_Len+1);

    float Diago, gapInA, gapInB;
    for(int i=1; i<=ChromA_Len; i++ ){
        for(int j=1; j<=ChromB_Len; j++){
            Diago = M[i-1][j-1] + s[i-1][j-1];
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
    printMatrix(*M, ChromA_Len+1, ChromB_Len+1);
    printMatrix((char *)Traceback, ChromA_Len+1, ChromB_Len+1);

    vector<int> indexA_aligned;
    vector<int> indexB_aligned;
    vector<float> score;
    char TracebackPointer;
    ROW_IDX = ChromA_Len;
    COL_IDX = ChromB_Len;

    bool OverlapAlignment = true;

    if(OverlapAlignment == true){
        // Overlap Alignment
        getOlapAlignStartIndices(*M, ChromA_Len+1, ChromB_Len+1, ROW_IDX, COL_IDX);
        //printMatrix(*M, ChromA_Len+1, ChromB_Len+1);
        //cout << ROW_IDX << " " << COL_IDX << endl;

        if(ROW_IDX != ChromA_Len){
            for (int i = ChromA_Len; i>ROW_IDX; i--){
                indexA_aligned.insert(indexA_aligned.begin(), i);
                indexB_aligned.insert(indexB_aligned.begin(), NA);
                score.insert(score.begin(), M[i][COL_IDX]);
            }
        } else if (COL_IDX != ChromB_Len){
            for (int j = ChromB_Len; j>COL_IDX; j--){
                indexA_aligned.insert(indexA_aligned.begin(), NA);
                indexB_aligned.insert(indexB_aligned.begin(), j);
                score.insert(score.begin(), M[ROW_IDX][j]);
            }
        }

        cout << endl;
        for (vector<float>::iterator it = score.begin(); it != score.end(); it++)
            cout << *it << " ";

        cout << endl;
        for (vector<int>::iterator it = indexA_aligned.begin(); it != indexA_aligned.end(); it++)
            cout << *it << " ";

        cout << endl;
        for (vector<int>::iterator it = indexB_aligned.begin(); it != indexB_aligned.end(); it++)
            cout << *it << " ";
        cout << endl;
    }

    TracebackPointer = Traceback[ROW_IDX][COL_IDX];
    cout << "Perform alignment" << endl;
    while(TracebackPointer != 'S'){
        // D: Diagonal, T: Top, L: Left
        if(TracebackPointer == 'D'){
            indexA_aligned.insert(indexA_aligned.begin(), ROW_IDX);
            indexB_aligned.insert(indexB_aligned.begin(), COL_IDX);
            score.insert(score.begin(), M[ROW_IDX][COL_IDX]);
            ROW_IDX = ROW_IDX-1;
            COL_IDX = COL_IDX-1;
        } else if(TracebackPointer == 'T') {
            indexA_aligned.insert(indexA_aligned.begin(), ROW_IDX);
            indexB_aligned.insert(indexB_aligned.begin(), NA);
            score.insert(score.begin(), M[ROW_IDX][COL_IDX]);
            ROW_IDX = ROW_IDX-1;
        } else {
            indexA_aligned.insert(indexA_aligned.begin(), NA);
            indexB_aligned.insert(indexB_aligned.begin(), COL_IDX);
            score.insert(score.begin(), M[ROW_IDX][COL_IDX]);
            COL_IDX = COL_IDX-1;
        }
        TracebackPointer = Traceback[ROW_IDX][COL_IDX];
        }
    cout << endl;
    for (vector<float>::iterator it = score.begin(); it != score.end(); it++)
        cout << *it << " ";

    cout << endl;
    for (vector<int>::iterator it = indexA_aligned.begin(); it != indexA_aligned.end(); it++)
        cout << *it << " ";

    cout << endl;
    for (vector<int>::iterator it = indexB_aligned.begin(); it != indexB_aligned.end(); it++)
        cout << *it << " ";
    cout << endl;

    return 0;
}
