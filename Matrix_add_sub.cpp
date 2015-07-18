#include<iostream>
#include<fstream>
#include<vector>
#include<cstdlib>

// Created by tarmizi adam, 18/7/2015. Simple snippet for two matrix addition and subtraction.
// Given two matrix A and B, the snippet produces C with the same dimension as A and B.

using namespace std;

vector< vector<double> > addMatrix(vector< vector<double> > &A, vector< vector<double> > &B);
vector< vector<double> > subMatrix(vector< vector<double> > &A, vector< vector<double> > &B);

int main()
{
    vector< vector<double> > A = {{1,1,2,3,5,8,13,21},
                                  {1,2,3,5,2,1,8,1},
                                  {1,1,1,1,1,1,1,1}};

    vector< vector<double> > B = {{1,1,1,1,1,1,1,1},
                                  {1,1,1,1,1,1,1,1},
                                  {1,1,1,1,1,1,1,1}};

    vector< vector<double> > C  = subMatrix(A,B);

    cout << C[0].size() << endl;

    cout << "A - B = " << endl;

    for(int i =0; i < C.size(); i++)
    {
        for(int j =0; j < C[0].size();j++)
        {
            cout << C[i][j] << ",";
        }

        cout << endl;
    }

    return 0;
}

vector< vector<double> > addMatrix(vector< vector<double> > &A, vector< vector<double> > &B)
{
    vector< vector<double> > C(A.size(), vector<double>(A[0].size(),0.0));

    if(A.size() != B.size() || A[0].size()!= B[0].size())
    {
        cout << "Error, Matrix dimension are not the same" << endl;
        exit(EXIT_FAILURE);
    }

    for(int i =0; i < A.size(); i++)
    {
        for(int j =0; j < A[0].size(); j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }

    return C;
}

vector< vector<double> > subMatrix(vector< vector<double> > &A, vector< vector<double> > &B)
{
     vector< vector<double> > C(A.size(), vector<double>(A[0].size(),0.0));

    if(A.size() != B.size() || A[0].size()!= B[0].size())
    {
        cout << "Error, Matrix dimension are not the same" << endl;
        exit(EXIT_FAILURE);
    }

    for(int i =0; i < A.size(); i++)
    {
        for(int j =0; j < A[0].size(); j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }

    return C;

}



