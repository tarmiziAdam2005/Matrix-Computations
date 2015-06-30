#include<iostream>
#include<vector>

// Created by Tarmizi Adam 30/6/2015. This simple snippet demostrates
// matrix vector (A x b) multiplication.

using namespace std;

void matrixVectMuliply(vector<double>&A, vector< vector<double> > &B);

int main()
{
    vector<double> X = {1,2,3,4};

    vector< vector<double> > Y = {{1,0,0,0},
                                    {2,1,0,0},
                                    {2,2,1,0},
                                    {1,2,2,1}};

    matrixVectMuliply(X, Y); // Multiply vector X with matrix Y.

    cout << endl;
    return 0;
}

void matrixVectMuliply(vector<double>&A, vector< vector<double> > &B)
{
    vector<double> Y(A.size(),0.0); // We must initialize our resultant matrix to 0
    int n = A.size();
    int m = A.size();
    int p = B[0].size();

    // Here we do some checking.
    // If the Column of A is not the same as the Row of B
    // multiplication cannot be done !
    if(A.size() != B.size())
    {
        cout << "Matrix dimension are not the same !" << endl;
        return;
    }

    // Matrix vector multiplication
    cout << "A x b = " << endl;

    for(int j = 0; j < p; j++)
    {
        double sum = 0;

        for(int k = 0; k < m; k++)
        {
            sum = sum + B[j][k]*A[k];
        }

        Y[j] = sum;
    }

    for(int i = 0; i < Y.size(); i++)
    {
        cout << Y[i] << endl;
    }
}
