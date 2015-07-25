#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>

//Created on 16/7/2015 by Tarmizi Adam. This simple snippet creates
// Update 25/7/2015
//Create an M-by-M sparse banded matrix by taking the columns of e (vector e) and placing them along the diagonals of D (vector of vector D).

//sparseDiag() create a main diagonal matrix with element e (vector e).

// This snippet can be used to create difference matrices (1st order, 2nd order etc) for signal processing.

// Dependencies: None
// Unlike the SparseDiagonalMatrix.cpp, this version does not depend on OpenCV Mat structure.
// The output is also slightly different.
// Pure C++

// I personally think this code is a little bit messy and not elegent. If you have a better way to
// generate the output, please, please, mail me at tarmizi_adam2005@yahoo.com

using namespace std;

template<class T> void printArray2D(vector< vector<T> > &I); // Function to show our created matrix
vector< vector<double> > bandedMat(vector<double> &e,int N);
vector< vector<double> > sparseDiag(vector<double> &Diag);
void absVal(vector<double> &i);

int main()
{
    int N = 3 ;
    vector<double> e = {-1,2,-1}; // first order difference
    //vector< vector<double> > B = bandedMat(e,N);

    // Proof ? see the banded matrix yourself !
    absVal(e);

    vector< vector<double> > F = sparseDiag(e);
    //vector< vector<double> > F =bandedMat(e,N);

    printArray2D(F);

    return 0;
}

template<class T> void printArray2D(vector< vector<T> > &I)
{
    // This is how we iterate using an iterator for 2d vectors
    typename vector< vector<T> >::iterator row; // Iterator for row of 2d vector
    typename vector<T>::iterator col; // Iterator for columns of 2d vector

    cout << "Matrix size: " << "[" << I.size() << "x" << I[0].size() << "]" << endl; // print the row and columns

    for(row = I.begin(); row !=I.end(); row++)
    {
        for(col = row->begin(); col != row->end(); col++)
        {
            cout << *col << ","; // print the value contained in each row of our 2d vector.
        }

        cout << endl;
    }
}

vector< vector<double> > bandedMat(vector<double> &e,int N)
{
    // Do some checking, Only tridiagonals are allowed.
    if(e.size()>3)
    {
        cout << "Only tridiagonals are allowed. Input vector to function must be only 3 elements." << endl;
        exit(EXIT_FAILURE);
    }

    vector< vector<double> > D(N-1, vector<double>(N-1,0.0));

    for(size_t i = 0; i < D.size(); i++)
    {
        for(size_t j = 0; j < D.size(); j++)
        {
            if(i == j)
            {
                D[i][j] = e[1]; // put value of e[1] as the main diagonal of D
            }
        }
    }

    int j =0; // index
    int k =0; // index

    for(size_t i =0; i < D.size()-1;i++)
    {
        D[i][j+1] = e[0]; // put value e[0] as the upper diagonal
        j = j+1;
    }

    for(size_t i =0; i < D.size()-1;i++)
    {
        D[i+1][k] = e[2]; // put value of e[2] as the lower diagonal
        k = k+1;
    }

    return D;
}

vector< vector<double> > sparseDiag(vector<double> &Diag)
{
    vector< vector<double> > D(Diag.size(), vector<double>(Diag.size(),0.0));

    for(size_t i = 0; i < D.size(); i++)
    {
        for(size_t j = 0; j < D.size(); j++)
        {
            if(i == j)
            {
                D[i][j] = Diag[j]; // put value of e[1] as the main diagonal of D
            }
        }
    }

        return D;
}

void absVal(vector<double> &sig)
{
    for(size_t i =0; i < sig.size(); i++)
    {
        sig[i] = std::abs(sig[i]);
    }
}
