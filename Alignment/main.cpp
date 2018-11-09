#include <iostream>
#include <stdio.h>
#include <vector>
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

template< class T>
void printMatrix(T *s, int ROW_SIZE, int COL_SIZE){
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            cout << *((s+i*COL_SIZE) + j) << " ";
        }
        cout << endl;
    }
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
        int indexA_aligned[5];
        int indexB_aligned[5];
        float score[5];
    };

    struct NeedleAlignObj
    {
        AlignedIndices Traceback[1];
        int M[ROW_SIZE][COL_SIZE];
        float Gap;
    };

    struct OLapAlignObj
    {
        AlignedIndices Traceback[1];
        int M[ROW_SIZE][COL_SIZE];
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
    getseqSimMat(seq1, seq2, Match, MisMatch, &s[0][0]); // getseqSimMat(seq1, seq2, Match, MisMatch, (float *)s);
    cout << "Similarity matrix is : " << endl;
    printMatrix(*s, seq1Len, seq2Len);
    cout << sizeof(s) / sizeof(s[0]) << endl;
    cout << sizeof(s[0]) / sizeof(s[0][0]) << endl;


    int ChromA_Len = sizeof(s) / sizeof(s[0]);
    int ChromB_Len = sizeof(s[0]) / sizeof(s[0][0]);
    float M[ChromA_Len+1][ChromB_Len+1] = {0};

    // enum TbPointer{STOP='S', T='T', D='D', L='L'};
    char Traceback[ChromA_Len+1][ChromB_Len+1] = {'S'};
    cout <<"Strta" <<endl;
    printMatrix((int *)Traceback, ChromA_Len+1, ChromB_Len+1);// How to change it?
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
    printMatrix(*M, ChromA_Len+1, ChromB_Len+1);
    cout <<"Strta" <<endl;
    printMatrix((char *)Traceback, ChromA_Len+1, ChromB_Len+1);
    printMatrix((int *)Traceback, ChromA_Len+1, ChromB_Len+1);

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



    /***
for(j in 2:(ChromB_Len+1)){
    for(i in 2:(ChromA_Len+1)){
      Diago <- M[i-1,j-1] + s[i-1, j-1]
      gapInA <- M[i-1,j] - gap
      gapInB <- M[i,j-1] - gap
      M[i,j] <- max(Diago, gapInA, gapInB)
      if(M[i,j] == Diago) Traceback[["TrM"]][i, j] <- "D" # D: Diagonal
      else if(M[i,j] == gapInA) Traceback[["TrM"]][i, j] <- "L" # L: Left
      else Traceback[["TrM"]][i, j] <- "T" # T: Top
    }}
    ***/
    int i;
    cout << "Hello !" << endl;
    int jj;
    cout << Match << " "<< MisMatch << " " << go << " " << ge << endl;
    cout << seq1 << " " << seq2 << endl;
    cout << 14%3 << endl;
    return 0;
}

