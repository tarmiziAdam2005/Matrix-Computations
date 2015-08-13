#include<iostream>
#include<vector>

// Created by tramizi adam on 13/8/2015. This simple program does a kronecker tensor product between
// two matrix A and B. The function is an imitation (hopefully) of MATLABs kron() function.

// Notes: The kronecker tensor product of two matrices A and B with sizes m x n and p x q produces a
// larger matrix C with size m*p x n*q.

// Dependencies: none

using namespace std;

vector< vector<double> > kron(vector< vector<double> > &A, vector< vector<double> > &B);

int main()
{
    vector< vector<double> > a = {{1.0,0.0,0.0,0.0},
                                  {0.0,1.0,0.0,0.0},
                                  {0.0,0.0,1.0,0.0},
                                  {0.0,0.0,0.0,1.0}}; // matrices a and b

    vector< vector<double> > x = {{1.0,-1.0},{-1.0,1.0}};

    auto T = kron(a,x);

    // Display the resulting matrix T
    for(size_t i  = 0; i < T.size(); i++)
    {
        for(size_t j = 0; j < T[0].size(); j++)
        {
            cout << T[i][j] << ",";
        }

        cout << endl;
    }

    return 0;
}

vector< vector<double> > kron(vector< vector<double> > &A, vector< vector<double> > &B)
{
    size_t krnProdRow = A.size()*B.size();
    size_t krnProdCol = A[0].size()*B[0].size();
    size_t nRowA = A.size();
    size_t nColA = A[0].size();
    size_t nRowB = B.size();
    size_t nColB = B[0].size();

    vector< vector<double> > krnProd(krnProdRow,vector<double>(krnProdCol,0.0));

    size_t i,j,k,l;
    double valA =0.0;
    double valB = 0.0;

    // Looping for kronecker tensor products of A and B
    // Trace the index of the loops to understand the kronecker tensor product
    for(i = 0; i <nRowA ; i++)
    {
        for(j = 0; j < nColA; j++)
        {
            valA = A[i][j];

            for(k = 0; k < nRowB; k++)
            {
                for(l = 0; l < nColB; l++)
                {
                    valB = B[k][l];
                    krnProd[nRowB*i+k][nColB*j+l] = valA * valB;
                }
            }
        }
    }

    return krnProd;
}

