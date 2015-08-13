#include<iostream>
#include<vector>
#include<cstdlib>

// Created by tramizi adam on 13/8/2015. This simple program verifies the kronecker tensor product.
// 1) do kronecker tensor product A and B producing C
// 2) multiply C with vector b

// check the results using another software to verify this program.

// Dependencies: none

using namespace std;

vector< vector<double> > kron(vector< vector<double> > &A, vector< vector<double> > &B);
void bandedMatVecMult(vector< vector<double> > &A, vector<double> &b);

int main()
{
    vector< vector<double> > a = {{1.0,0.0,0.0,0.0},
                                  {0.0,1.0,0.0,0.0},
                                  {0.0,0.0,1.0,0.0},
                                  {0.0,0.0,0.0,1.0}}; // matrices a and b

    vector< vector<double> > x = {{1.0,-1.0},{-1.0,1.0}};

    vector<double> b = {2.0,5.0,8.0,10.0,1.0,2.0,3.0,4.0};

    vector< vector<double> > T = kron(a,x);

    bandedMatVecMult(T,b);

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

void bandedMatVecMult(vector< vector<double> > &A, vector<double> &b)
{
    vector<double>diags(A.size(),0.0);
    vector<double>upperBand(A.size(),0.0);
    vector<double>lowerBand(A.size(),0.0);

    vector<double> resUpperBand(upperBand.size(),0.0);
    vector<double> reslowerBand(lowerBand.size(),0.0);
    vector<double> resDiags(A.size(),0.0);
    vector<double> result(A.size(),0.0);

    size_t rows = A.size();
    size_t cols = A[0].size();
    size_t j =0;

    if(A[0].size() != b.size())
    {
        cout << "Multiplication of bad sizes...";
        exit(EXIT_FAILURE);
    }


    for(size_t i =0; i < rows; i++)
    {
        if(i == j)
        {
            diags[i] = A[i][j]; // Save the diagonals of A only....
        }
        j++;
    }

    //cout << "Upper band:";
    for(size_t i =0; i < upperBand.size()-1; i++)
    {
        upperBand[i] = A[i][i+1]; // Save the diagonal upper band of A only....
    }

    //cout << "Lower band:";
    for(size_t i =1; i < lowerBand.size(); i++)
    {
        lowerBand[i] = A[i][i-1]; // Save the diagonals lower band of A only....
    }

    //Upper band multiplication with b

    for(size_t i = 0; i < upperBand.size()-1; i++)
    {
        resUpperBand[i] = upperBand[i]*b[i+1];
    }

    //lower band multiplication with b
    reslowerBand[0] =0.0;
    for(size_t i = 1; i < lowerBand.size(); i++)
    {
        //cout << lowerBand[i] << endl;
        reslowerBand[i] = lowerBand[i]*b[i-1];
    }

    //Diagonal multiplication with b
    for(size_t i =0; i < diags.size(); i++)
    {
        resDiags[i] = diags[i]*b[i];
    }

    // Result of multiplication vector. Sum of the tri diagonals above

    for(size_t i = 0; i < result.size(); i++)
    {
        result[i] = resDiags[i]+reslowerBand[i]+resUpperBand[i];
    }

    cout <<"result: " << endl;
    for(size_t i =0; i < result.size(); i++)
    {
        cout << result[i] << endl;
    }

}
