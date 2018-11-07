#include <iostream>
using namespace std;


void newf(int &sim[N][M]){
    sim[N-1][M-1]=0;
}

void getseqSimMat(string seq1, string seq2, float Match, float MisMatch, float (&s)[rows][cols]){
    ROW_SIZE = seq1.size();
    COL_SIZE = seq2.size();
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            seq1[i] == seq2[j] ? s[i][j] = Match : s[i][j] = MisMatch;
        }
    }
}

void printLengthOfSeq(string seq){
    int ROW_SIZE;
    ROW_SIZE = seq.size();
    cout << "The length of sequence is " << seq[2] << endl;
}

int main()
{
    int ROW_SIZE = 3; // If I do not do const, it doesn't pass on to struct and throws error
    int COL_SIZE = 4; // that array bound is not an integer constant before ']'.

    int m[ROW_SIZE][COL_SIZE]=
        {
            {1,2,3,4},
            {5,6,7,8},
            {13,14,15,16}
        };


    int VecMat[ROW_SIZE*COL_SIZE]= {1,5,13, 2,6,14, 3,7,15, 4,8,16};

    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            cout << m[i][j] << " ";
        }
        cout << endl;
    }
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
    /***
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

    ***/

    //NeedleAlignObj doGlobalAlignment(float s, float gap){
    //}
    float Match=10, MisMatch=-2, go=22, ge=7;
    string seq1 = "GCAT";
    string seq2 = "CAGTG";
    float s[4][5];
    printLengthOfSeq(seq1);
    printLengthOfSeq(seq2);
    //getseqSimMat(seq1, seq2, Match, MisMatch, s);


    int i;
    cout << "Hello !" << endl;
    int jj;
    cout << Match << " "<< MisMatch << " " << go << " " << ge << endl;
    cout << seq1 << " " << seq2 << endl;
    cout << 14%3 << endl;
    return 0;
}

