#include <iostream>
#include <vector>

// Created by Tarmizi Adam 22/6/2015. A simple program that initializes a 2d vector (matrix)
// and prints the matrix. The function printArray2D is a template function thus, it can print out
// any type of matrix created (float, doubles, int)

// Input: Matrix to be printed
// Output: Displays the matrix to the user

// Dependencies: None

using namespace std;

template<class T> void printArray2D(vector< vector<T> > &I); // function template to print our matrix. Any type can be used
//void printArray2D(vector< vector<mpz_int>> &I)

int main()
{

   int N = 4;
   vector< vector<float> > A(N, vector<float>(N,0)); // Initialize A to zeros (example)

   printArray2D(A);

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
