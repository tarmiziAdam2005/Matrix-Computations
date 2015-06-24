#include<iostream>
#include<vector>

// Created by Tarmizi Adam 24/6/2015. This simple snippet demostrates
// matrix multiplication.

using namespace std;

void matrixMuliply(vector< vector<double> > &A, vector< vector<double> > &B);

int main()
{
    vector< vector<double> > X = {{4,2,3},
                                    {2,2,3},
                                    {1,3,10},
                                    {3,6,6}};

    vector< vector<double> > Y = {{1,0,0,0},
                                    {2,1,0,0},
                                    {2,2,1,0},
                                    {1,2,2,1}};

    matrixMuliply(Y, Y); // Multiply B and A.

    cout << endl;
    return 0;
}

void matrixMuliply(vector< vector<double> > &A, vector< vector<double> > &B)
{
    vector< vector<double> > Y(A.size(),vector<double>(B[0].size(),0.0)); // We must initialize our resultant matrix to 0
    int n = A.size();
    int m = A[0].size();
    int p = B[0].size();

    // Here we do some checking.
    // If the Column of A is not the same as the Row of B
    // multiplication cannot be done !
    if(A[0].size() != B.size())
    {
        cout << "Matrix dimension are not the same !" << endl;
        return;
    }

    // Matrix multiplication
    cout << "A x B = " << endl;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < p; j++)
        {
            double sum = 0;

            for(int k = 0; k < m; k++)
            {
                sum = sum + A[i][k]*B[k][j];
            }

            Y[i][j] = sum;
            cout << Y[i][j] <<",";
        }
        cout << endl;
    }
}
