#include<iostream>
#include<vector>
#include<cstdlib>

//Created on 18/7/2015 by Tarmizi Adam. This simple snippet creates
//an M-by-M sparse matrix by taking the elements of e (vector e) and placing them along the diagonals of D (vector of vector D).

// Dependencies: None

// comments bugs or questions, mail me at tarmizi_adam2005@yahoo.com

using namespace std;

template<class T> void printArray2D(vector< vector<T> > &I); // Function to show our created matrix
vector< vector<double> > sparseDiag(vector<double> &Diag);

int main()
{
    vector<double> e = {-1,2,-9,7,6,-3,8,-4};

    vector< vector<double> > F = sparseDiag(e);

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

vector< vector<double> > sparseDiag(vector<double> &Diag)
{
    vector< vector<double> > D(Diag.size(), vector<double>(Diag.size(),0.0));

    for(size_t i = 0; i < D.size(); i++)
    {
        for(size_t j = 0; j < D.size(); j++)
        {
            if(i == j)
            {
                D[i][j] = Diag[j]; // put values of Diag[1] as the main diagonal of D
            }
        }
    }

        return D;
}
