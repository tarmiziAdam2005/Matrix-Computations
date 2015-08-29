#include<iostream>
#include<vector>

// Created on 30/8/2015, Tarmizi Adam. Resize 2D matrix to 1D vector.
// Similar to doing (:) in MATLAB. i.e. a = a(:).

using namespace std;

int main()
{
    vector< vector<double> > a(3,vector<double>(3,0.0));
    int N = a.size()*a.size();

    a = {{1,2,3},
         {4,5,5},
         {6,6,6}};

    int n = a[0].size();

    vector<double> resized(N,0.0);

    for(int i = 0 ; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            resized[j*n + i] = a[i][j]; // resized to [j*n + i] sized vector.
        }

    }

    for(int i = 0; i < resized.size(); i++)
    {
        cout << resized[i] << "," << endl;
    }

    return 0;
}
