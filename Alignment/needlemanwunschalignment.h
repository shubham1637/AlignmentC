#ifndef NEEDLEMANWUNSCHALIGNMENT
#define NEEDLEMANWUNSCHALIGNMENT

#endif // NEEDLEMANWUNSCHALIGNMENT

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
TracebackPointer = Traceback[ChromA_Len][ChromB_Len];
ROW_IDX = ChromA_Len;
COL_IDX = ChromB_Len;
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
        indexB_aligned.insert(indexB_aligned.begin(), -1);
        score.insert(score.begin(), M[ROW_IDX][COL_IDX]);
        ROW_IDX = ROW_IDX-1;
    } else {
        indexA_aligned.insert(indexA_aligned.begin(), -1);
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
