#include<iostream>
#include<vector>

//Created on 16/7/2015 by Tarmizi Adam. This simple snippet creates
//an m-by-n sparse matrix by taking the columns of e (vector e) and placing them along the diagonals of D (vector of vector D).

// This snippet can be used to create difference matrices (1st order, 2nd order etc) for signal processing.

// Dependencies: None
// Unlike the SparseDiagonalMatrix.cpp, this version does not depend on OpenCV Mat structure.
// Pure C++

using namespace std;

template<class T> void printArray2D(vector< vector<T> > &I); // Function to show our created matrix

int main()
{
    int N = 15 ;
    vector< vector<double> > D(N-2, vector<double>(N,0.0));
    double val = 0.0;

    vector<double> e = {1,-2,1}; // first order difference

     for(int i = 0; i < D.size(); i++)
    {
        for(int j = 0; j < e.size(); j++)
        {
            val = e[j];
            D[i][i+j] = val; // Put along the diagonal of matrix D. Note the index.
        }
    }

    // Proof ? see the banded matrix yourself !
    printArray2D(D);

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
