#include <iostream>
using namespace std;

int main()
{
    int ROW_SIZE = 3;
    int COL_SIZE = 4;

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


    int i;
    cout << "Hello !" << endl;
    int jj;
    cout << 14/3 << endl;
    cout << 14%3 << endl;
    return 0;
}

