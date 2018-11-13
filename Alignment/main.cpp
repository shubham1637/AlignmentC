#include <iostream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <limits>
#define NA -1

//using namespace std;

enum TracebackType {SS = 0, DM = 1, DA = 2, DB = 3, TM = 4, TA = 5, TB = 6, LM = 7, LA = 8, LB = 9};
std::ostream& operator<<(std::ostream& out, const TracebackType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(SS);
        PROCESS_VAL(DM);
        PROCESS_VAL(DA);
        PROCESS_VAL(DB);
        PROCESS_VAL(TM);
        PROCESS_VAL(TA);
        PROCESS_VAL(TB);
        PROCESS_VAL(LM);
        PROCESS_VAL(LA);
        PROCESS_VAL(LB);
    }
#undef PROCESS_VAL
    return out << s;
}
enum tbJump {M = 0, A = 1, B = 2};
struct AffineAlignObj
{
    float* M;
    float* A;
    float* B;
    TracebackType* Traceback;
    int signalA_len; // stack allocation
    int signalB_len;
    float GapOpen;
    float GapExten;
    bool FreeEndGaps;

    // Not a default constructor
    AffineAlignObj(int ROW_SIZE, int COL_SIZE)
    {
        M = new float[ROW_SIZE * COL_SIZE]; // heap allocation
        A = new float[ROW_SIZE * COL_SIZE]; // heap allocation
        B = new float[ROW_SIZE * COL_SIZE]; // heap allocation
        Traceback = new TracebackType[3*ROW_SIZE * COL_SIZE]; // heap allocation
        signalA_len = ROW_SIZE-1;
        signalB_len = COL_SIZE-1;
    }

    // Rule 1 Copy constructor
    AffineAlignObj(const AffineAlignObj &other)
       {
          signalA_len = other.signalA_len;
          signalB_len = other.signalB_len;
          GapOpen = other.GapOpen;
          GapExten = other.GapExten;
          FreeEndGaps = other.FreeEndGaps;
          M = new float[(signalA_len+1)*(signalB_len+1)];
          std::memcpy(M, other.M, sizeof(float) * (signalA_len+1)*(signalB_len+1));
          A = new float[(signalA_len+1)*(signalB_len+1)];
          std::memcpy(A, other.A, sizeof(float) * (signalA_len+1)*(signalB_len+1));
          B = new float[(signalA_len+1)*(signalB_len+1)];
          std::memcpy(B, other.B, sizeof(float) * (signalA_len+1)*(signalB_len+1));
          Traceback = new TracebackType[3*(signalA_len+1)*(signalB_len+1)];
          std::memcpy(Traceback, other.Traceback, sizeof(TracebackType) * 3 * (signalA_len+1)*(signalB_len+1));
       }

    // Rule 2 Copy assignment operator
    AffineAlignObj& operator=(const AffineAlignObj& other)
    {
        if(this == &other) return *this; // handling of self assignment.
        delete[] M; // freeing previously used memory
        delete[] A;
        delete[] B;
        delete[] Traceback;
        signalA_len = other.signalA_len;
        signalB_len = other.signalB_len;
        GapOpen = other.GapOpen;
        GapExten = other.GapExten;
        FreeEndGaps = other.FreeEndGaps;
        M = new float[(signalA_len+1)*(signalB_len+1)];
        std::memcpy(M, other.M, sizeof(float) * (signalA_len+1)*(signalB_len+1));
        A = new float[(signalA_len+1)*(signalB_len+1)];
        std::memcpy(A, other.A, sizeof(float) * (signalA_len+1)*(signalB_len+1));
        B = new float[(signalA_len+1)*(signalB_len+1)];
        std::memcpy(B, other.B, sizeof(float) * (signalA_len+1)*(signalB_len+1));
        Traceback = new TracebackType[3*(signalA_len+1)*(signalB_len+1)];
        std::memcpy(Traceback, other.Traceback, sizeof(TracebackType) * 3 * (signalA_len+1)*(signalB_len+1));
        return *this;
    }

    // Rule 3 Not a default destructor
    ~AffineAlignObj()
    {
        delete[] M; // since we declared with new, manually clear memory from heap.
        delete[] A;
        delete[] B;
        delete[] Traceback;
    }
};

void getseqSimMat(std::string seq1, std::string seq2, float Match, float MisMatch, float *s){
    int ROW_SIZE = seq1.size();
    int COL_SIZE = seq2.size();
    for(int i = 0; i < ROW_SIZE; i++){
        for(int j = 0; j < COL_SIZE; j++){
            seq1[i] == seq2[j] ? *((s+i*COL_SIZE) + j) = Match : *((s+i*COL_SIZE) + j) = MisMatch;
        }
    }
}

void printLengthOfSeq(std::string seq){
    std::cout << "The length of sequence is " << seq.size() << std::endl;
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

template<class T>
T getOlapAffineAlignStartIndices(T *MatrixM, T *MatrixA, T *MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName){
    T affineAlignmentScore;
    float maxScore = -std::numeric_limits<float>::infinity();
    int MaxRowIndex, MaxColIndex;
    for(int i = 0; i < ROW_SIZE; i++){
        if(*(MatrixM+i*COL_SIZE+COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            MatrixName = M;
            maxScore = *(MatrixM+i*COL_SIZE+COL_SIZE-1);
        } else if(*(MatrixA+i*COL_SIZE+COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            MatrixName = A;
            maxScore = *(MatrixA+i*COL_SIZE+COL_SIZE-1);
        } else if(*(MatrixB+i*COL_SIZE+COL_SIZE-1) >= maxScore){
            MaxRowIndex = i;
            MaxColIndex = COL_SIZE-1;
            MatrixName = B;
            maxScore = *(MatrixB+i*COL_SIZE+COL_SIZE-1);
        }
    }
    for (int j = 0; j < COL_SIZE; j++){
        if(*(MatrixM+(ROW_SIZE-1)*COL_SIZE+j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            MatrixName = M;
            maxScore = *(MatrixM+(ROW_SIZE-1)*COL_SIZE+j);
        } else if(*(MatrixA+(ROW_SIZE-1)*COL_SIZE+j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            MatrixName = A;
            maxScore = *(MatrixA+(ROW_SIZE-1)*COL_SIZE+j);
        } else if(*(MatrixB+(ROW_SIZE-1)*COL_SIZE+j) >= maxScore){
            MaxRowIndex = ROW_SIZE-1;
            MaxColIndex = j;
            MatrixName = B;
            maxScore = *(MatrixB+(ROW_SIZE-1)*COL_SIZE+j);
        }
    }
    OlapStartRow = MaxRowIndex;
    OlapStartCol = MaxColIndex;
    affineAlignmentScore = maxScore;
    return affineAlignmentScore;
}

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

void doAffineAlignment(float *s, int signalA_len, int signalB_len, float go, float ge, bool OverlapAlignment, AffineAlignObj &affineAlignObj){
    affineAlignObj.FreeEndGaps = OverlapAlignment;
    affineAlignObj.GapOpen = go;
    affineAlignObj.GapExten = ge;

    // Initialize first row and first column for global and overlap alignment.
    float Inf = std::numeric_limits<float>::infinity();
    for(int i = 0; i<=signalA_len; i++){
        *(affineAlignObj.M+i*(signalB_len+1)+0) = -Inf;
        *(affineAlignObj.B+i*(signalB_len+1)+0) = -Inf;
        *(affineAlignObj.Traceback + 0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0) = SS; //STOP
        *(affineAlignObj.Traceback + 2*((signalA_len+1)*(signalB_len+1)) + i*(signalB_len+1)+0) = SS; //STOP
    }
    for(int j = 0; j<=signalB_len; j++){
        *(affineAlignObj.M+0*(signalB_len+1)+j) = -Inf;
        *(affineAlignObj.A+0*(signalB_len+1)+j) = -Inf;
        *(affineAlignObj.Traceback + 0*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j) = SS; //STOP
        *(affineAlignObj.Traceback + 1*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j) = SS; //STOP
    }
    *(affineAlignObj.M+0*(signalB_len+1)+0) = 0;
    if(affineAlignObj.FreeEndGaps == true){
        for(int i = 1; i<=signalA_len; i++){
            *(affineAlignObj.A+i*(signalB_len+1)+0) = 0;
            *(affineAlignObj.Traceback + 1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0) = TA; //TOP A
        }
        for(int j = 1; j<=signalB_len; j++){
            *(affineAlignObj.B+0*(signalB_len+1)+j) = 0;
            *(affineAlignObj.Traceback + 2*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j) = LB; //LEFT B
        }
    } else {
        for(int i = 1; i<=signalA_len; i++){
            *(affineAlignObj.A+i*(signalB_len+1)+0) = -(i-1)*ge - go;
            *(affineAlignObj.Traceback + 1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+0) = TA; //TOP A
        }
        for(int j = 1; j<=signalB_len; j++){
            *(affineAlignObj.B+0*(signalB_len+1)+j) = -(j-1)*ge - go;
            *(affineAlignObj.Traceback + 2*((signalA_len+1)*(signalB_len+1))+ 0*(signalB_len+1)+j) = LB; //LEFT B
        }
    }

    // Perform dynamic programming for affine alignment
    float Diago, gapInA, gapInB;
    for(int i=1; i<=signalA_len; i++){
        for(int j=1; j<=signalB_len; j++){
            float sI_1J_1 = *((s+(i-1)*signalB_len) + j-1);
            Diago = *(affineAlignObj.M+(i-1)*(signalB_len+1)+j-1) + sI_1J_1;
            gapInA = *(affineAlignObj.A+(i-1)*(signalB_len+1)+j-1) + sI_1J_1;
            gapInB = *(affineAlignObj.B+(i-1)*(signalB_len+1)+j-1) + sI_1J_1;

            // Calculate recursively for matched alignment
            if(Diago>=gapInA && Diago>=gapInB){
                *(affineAlignObj.Traceback + 0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = DM; // DM: Diagonal TrM
                *(affineAlignObj.M+i*(signalB_len+1)+j) = Diago;
            }
            else if (gapInA>=Diago && gapInA>=gapInB){
                *(affineAlignObj.Traceback + 0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = DA; // DA: Diagonal TrA
                *(affineAlignObj.M+i*(signalB_len+1)+j) = gapInA;
            }
            else{
                *(affineAlignObj.Traceback + 0*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = DB; // DB: Diagonal TrB
                *(affineAlignObj.M+i*(signalB_len+1)+j) = gapInB;
            }

            // Calculate recursively for gap in signalB
            float AfromM = *(affineAlignObj.M+(i-1)*(signalB_len+1)+j) - go;
            float AfromA = *(affineAlignObj.A+(i-1)*(signalB_len+1)+j) - ge;
            float AfromB = *(affineAlignObj.B+(i-1)*(signalB_len+1)+j) - go;
            if(AfromM >= AfromA && AfromM >= AfromB){
                *(affineAlignObj.Traceback + 1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = TM; // TM: Top TrM
                *(affineAlignObj.A+i*(signalB_len+1)+j) = AfromM;
            }
            else if (AfromA >= AfromM && AfromA >= AfromB){
                *(affineAlignObj.Traceback + 1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = TA; // TA: Top TrA
                *(affineAlignObj.A+i*(signalB_len+1)+j) = AfromA;
            }
            else{
                *(affineAlignObj.Traceback + 1*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = TB; // TB: Top TrB
                *(affineAlignObj.A+i*(signalB_len+1)+j) = AfromB;
            }

            // Calculate recursively for gap in signalA
            float BfromM = *(affineAlignObj.M+i*(signalB_len+1)+j-1) - go;
            float BfromA = *(affineAlignObj.A+i*(signalB_len+1)+j-1) - go;
            float BfromB = *(affineAlignObj.B+i*(signalB_len+1)+j-1) - ge;
            if(BfromM >= BfromA && BfromM >= BfromB){
                *(affineAlignObj.Traceback + 2*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = LM; // LM: Left TrM
                *(affineAlignObj.B+i*(signalB_len+1)+j) = BfromM;
            }
            else if (BfromA >= BfromM && BfromA >= BfromB){
                *(affineAlignObj.Traceback + 2*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = LA; // LA: Left TrA
                *(affineAlignObj.B+i*(signalB_len+1)+j) = BfromA;
            }
            else{
                *(affineAlignObj.Traceback + 2*((signalA_len+1)*(signalB_len+1))+ i*(signalB_len+1)+j) = LB; // LB: Left TrB
                *(affineAlignObj.B+i*(signalB_len+1)+j) = BfromB;
            }
        }
    }
}

AlignedIndices getAffineAlignedIndices(AffineAlignObj affineAlignObj){
    AlignedIndices alignedIdx;
    char TracebackPointer;
    tbJump MatName;
    float affineAlignmentScore;
    int ROW_IDX = affineAlignObj.signalA_len;
    int COL_IDX = affineAlignObj.signalB_len;
    int ROW_SIZE = (affineAlignObj.signalA_len)+1;
    int COL_SIZE = (affineAlignObj.signalB_len)+1;

    if(affineAlignObj.FreeEndGaps == true){
        // Overlap Alignment
        affineAlignmentScore = getOlapAffineAlignStartIndices((float *)affineAlignObj.M, (float *)affineAlignObj.A, (float *)affineAlignObj.B, ROW_SIZE, COL_SIZE, ROW_IDX, COL_IDX, MatName);
        if(ROW_IDX != affineAlignObj.signalA_len){
            for (int i = affineAlignObj.signalA_len; i>ROW_IDX; i--){
                alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), i);
                alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
                alignedIdx.score.insert(alignedIdx.score.begin(), affineAlignmentScore);
            }
        } else if (COL_IDX != affineAlignObj.signalB_len){
            for (int j = affineAlignObj.signalB_len; j>COL_IDX; j--){
                alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
                alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), j);
                alignedIdx.score.insert(alignedIdx.score.begin(), affineAlignmentScore);
            }
        }
    } else {
        // Global Alignment
        float Mscore = *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX);
        float Ascore = *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX);
        float Bscore = *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX);
        if (Mscore >= Ascore && Mscore >= Bscore){
            affineAlignmentScore = Mscore;
            MatName = M;
        } else if(Ascore >= Mscore && Ascore>= Bscore) {
            affineAlignmentScore = Ascore;
            MatName = A;
        } else{
            affineAlignmentScore = Bscore;
            MatName = B;
        }
    }

    alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
    alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
    alignedIdx.score.insert(alignedIdx.score.begin(), affineAlignmentScore);
    TracebackPointer = *(affineAlignObj.Traceback+MatName*ROW_SIZE*COL_SIZE+ROW_IDX*COL_SIZE+COL_IDX);
    // Traceback path and align row indices to column indices.
    while(TracebackPointer != SS){
        // D: Diagonal, T: Top, L: Left
        switch(TracebackPointer){
        // In the code below, we are appending future values. Because, once we go to the M,A or B matrix.
        // we will not be able to tell which matrix we are currently in.
        case DM:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = M;
            COL_IDX = COL_IDX-1;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(),  *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case DA:
        {
            ROW_IDX = ROW_IDX-1;
            COL_IDX = COL_IDX-1;
            MatName = A;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(),  *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case DB:
        {
            ROW_IDX = ROW_IDX-1;
            COL_IDX = COL_IDX-1;
            MatName = B;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(),  *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case TM:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = M;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case TA:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = A;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case TB:
        {
            ROW_IDX = ROW_IDX-1;
            MatName = B;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), ROW_IDX);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), NA);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case LM:
        {
            COL_IDX = COL_IDX-1;
            MatName = M;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.M+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case LA:
        {
            COL_IDX = COL_IDX-1;
            MatName = A;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.A+ROW_IDX*COL_SIZE+COL_IDX));
        }
        case LB:
        {
            COL_IDX = COL_IDX-1;
            MatName = B;
            alignedIdx.indexA_aligned.insert(alignedIdx.indexA_aligned.begin(), NA);
            alignedIdx.indexB_aligned.insert(alignedIdx.indexB_aligned.begin(), COL_IDX);
            alignedIdx.score.insert(alignedIdx.score.begin(), *(affineAlignObj.B+ROW_IDX*COL_SIZE+COL_IDX));
        }
        }
        TracebackPointer = *(affineAlignObj.Traceback+MatName*ROW_SIZE*COL_SIZE+ROW_IDX*COL_SIZE+COL_IDX);
        }

    return alignedIdx;
}


int main()
{
    float Match=10, MisMatch=-2, go=22, ge=7, gap=go;
    std::string seq1 = "GCAT";
    std::string seq2 = "CAGTG";
    int seq1Len = seq1.size();
    int seq2Len = seq2.size();
    float s[seq1Len][seq2Len];
    initializeMatrix(*s, 0, seq1Len, seq2Len);
    getseqSimMat(seq1, seq2, Match, MisMatch, &s[0][0]); // getseqSimMat(seq1, seq2, Match, MisMatch, (float *)s);
    std::cout << "Similarity matrix is : " << std::endl;
    printMatrix(*s, seq1Len, seq2Len);
    int signalA_len = sizeof(s) / sizeof(s[0]);
    int signalB_len = sizeof(s[0]) / sizeof(s[0][0]);
    bool OverlapAlignment = true;

    AlignObj alignObj;
    alignObj = doAlignment(*s, seq1Len, seq2Len, gap, OverlapAlignment);
    AffineAlignObj affineAlignObj(seq1Len+1, seq2Len+1); // What if this length is different than used inside the function.
    doAffineAlignment(*s, seq1Len, seq2Len, go, ge, OverlapAlignment, affineAlignObj);
    std::cout << "Out the loop " << std::endl;

//    AffineAlignObj affineAlignObj2(seq1Len+1, seq2Len+1);
//    affineAlignObj2 = affineAlignObj;
//    doAffineAlignment(*s, seq1Len, seq2Len, 3, 3, OverlapAlignment, affineAlignObj);

    std::cout << "M matrix is : " << std::endl;
    printMatrix((float *)affineAlignObj.M, signalA_len+1, signalB_len+1);
    printMatrix((float *)affineAlignObj.A, signalA_len+1, signalB_len+1);
    printMatrix((float *)affineAlignObj.B, signalA_len+1, signalB_len+1);
    printMatrix((TracebackType *)affineAlignObj.Traceback, signalA_len+1, signalB_len+1);
    printMatrix((TracebackType *)affineAlignObj.Traceback+1*((signalA_len+1)*(signalB_len+1)), signalA_len+1, signalB_len+1);
    printMatrix((TracebackType *)affineAlignObj.Traceback+2*((signalA_len+1)*(signalB_len+1)), signalA_len+1, signalB_len+1);

//    std::cout << "M matrix is : " << std::endl;
//    printMatrix((float *)affineAlignObj2.M, signalA_len+1, signalB_len+1);
//    printMatrix((float *)affineAlignObj2.A, signalA_len+1, signalB_len+1);
//    printMatrix((float *)affineAlignObj2.B, signalA_len+1, signalB_len+1);

    AlignedIndices alignedIdx;
    alignedIdx = getAffineAlignedIndices(affineAlignObj);

    std::cout << std::endl;
    for (std::vector<float>::iterator it = alignedIdx.score.begin(); it != alignedIdx.score.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;
    for (std::vector<int>::iterator it = alignedIdx.indexA_aligned.begin(); it != alignedIdx.indexA_aligned.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;
    for (std::vector<int>::iterator it = alignedIdx.indexB_aligned.begin(); it != alignedIdx.indexB_aligned.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;

    return 0;
}
