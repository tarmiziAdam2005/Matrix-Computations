#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>

//Created on 16/7/2015 by Tarmizi Adam. This simple snippet creates
// Update 25/7/2015
// Update26/7/2015
//Create an M-by-M sparse banded matrix by taking the columns of e (vector e) and placing them along the diagonals of D (vector of vector D).

//sparseDiag() create a main diagonal matrix with element e (vector e).

// This code does multiplication of sparse tri-diagonal matrix with a vector i.e. A*b.
// The zero entries in the sparse tri-diagonal matrix is discarded for computation efficiency.
// upper, lower and main diagonals are stored in three separate vectors. Then the computations are done separately.
// By this, zero entries of the sparse tri-diagonal matrice are discarded in the computation.
// For details, study the functions sparseDiagMatVecMult() and bandedMatVecMult();

// Dependencies: None
// Pure C++

// I personally think this code is a little bit messy and not elegent. If you have a better way to
// generate the output, please, please, mail me at tarmizi_adam2005@yahoo.com

using namespace std;

template<class T> void printArray2D(vector< vector<T> > &I); // Function to show our created matrix
vector< vector<double> > bandedMat(vector<double> &e,int N);
vector< vector<double> > sparseDiag(vector<double> &Diag);
vector<double> sparseDiagMatVecMult(vector< vector<double> > &A, vector<double> &b);
void bandedMatVecMult(vector< vector<double> > &A, vector<double> &b);



int main()
{
    int N = 5 ;
    vector<double> e = {1,2,1}; // test vector
    vector<double> tst = {2,5,8,10};// test vector

    //vector< vector<double> > B = bandedMat(e,N);
    //vector< vector<double> > F = sparseDiag(e);

    vector< vector<double> > F =bandedMat(e,N);

    //printArray2D(F);
    //cout << tst.size() << endl;
    //sparseDiagMatVecMult(F,e);
    //bandedMatVecMult(F,e);
    bandedMatVecMult(F,tst); // Sparse tri-diagonal multiply with vector (tst). A*b

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

vector<double> sparseDiagMatVecMult( vector< vector<double> > &A, vector<double> &b)
{
    vector<double>diags(A.size(),0.0);
    vector<double> res(diags.size(),0.0);

    int rows = A.size();
    int cols = A[0].size();
    int j =0;

    for(int i =0; i < rows; i++)
    {
        if(i == j)
        {
            diags[i] = A[i][j]; // Save the diagonals of A only....
        }
        j++;
    }

    if(diags.size()!=b.size())
    {
        cout << "Vector size not the same. Unable to do multiplication...";
        exit(EXIT_FAILURE);
    }

    for(int j =0; j < diags.size(); j++)
    {
        res[j] = diags[j]*b[j]; // Do the multiplication of diagonals and the vector.
    }

    // For display, comment out if not relevant
    /*for(int k =0; k < res.size(); k++)
    {
        cout << res[k] << endl;
    }*/

    return res;
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

    int rows = A.size();
    int cols = A[0].size();
    int j =0;

    for(int i =0; i < rows; i++)
    {
        if(i == j)
        {
            diags[i] = A[i][j]; // Save the diagonals of A only....
        }
        j++;
    }

    //cout << "Upper band:";
    for(int i =0; i < upperBand.size()-1; i++)
    {
        upperBand[i] = A[i][i+1]; // Save the diagonal upper band of A only....
    }

    //cout << "Lower band:";
    for(int i =1; i < lowerBand.size(); i++)
    {
        lowerBand[i] = A[i][i-1]; // Save the diagonals lower band of A only....
    }

    //Upper band multiplication with b

    for(int i = 0; i < upperBand.size()-1; i++)
    {
        resUpperBand[i] = upperBand[i]*b[i+1];
    }

    //lower band multiplication with b
    reslowerBand[0] =0.0;
    for(int i = 1; i < lowerBand.size(); i++)
    {
        //cout << lowerBand[i] << endl;
        reslowerBand[i] = lowerBand[i]*b[i-1];
    }

    //Diagonal multiplication with b
    for(int i =0; i < diags.size(); i++)
    {
        resDiags[i] = diags[i]*b[i];
    }

    // Result of multiplication vector. Sum of the tri diagonals above

    for(int i = 0; i < result.size(); i++)
    {
        result[i] = resDiags[i]+reslowerBand[i]+resUpperBand[i];
    }

    cout <<"lowerBand: " << endl;
    for(int i =0; i < lowerBand.size(); i++)
    {
        cout << lowerBand[i] << endl;
    }

    cout <<"result: " << endl;
    for(int i =0; i < result.size(); i++)
    {
        cout << result[i] << endl;
    }

}
