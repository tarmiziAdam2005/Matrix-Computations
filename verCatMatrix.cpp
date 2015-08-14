#include<iostream>
#include<vector>

//Created on 14/8/2015 by Tarmizi. This is a simple program to vertically concatenate matrices.
// Function arguments: res     --> the resulting matrix
//                     A       --> A matrix (2D vector) that you want to vertically concatenate

// See the example on how to use it

using namespace std;

template<class T> void verCatMatrix(vector< vector<T> > &res, vector< vector<T> > &A);

int main()
{
    vector< vector<int> > res;
    vector< vector<int> >a = {{4,2,3},
                              {2,2,2}};

    vector< vector<int> >b = {{4,5,6},
                              {6,6,6}};

    vector< vector<int> >t = {{8,8,8},
                              {10,10,10}};

    // vertically concatenate matrices a, b, and t
    verCatMatrix(res,a);
    verCatMatrix(res,b);
    verCatMatrix(res,t);

    // Display the resulting matrix
    for(int size_t = 0; i < res.size(); i++)
    {

        for(int size_t = 0; j < res[0].size(); j++)
        {
            cout << res[i][j] << ",";
        }

        cout << endl;
    }

    return 0;
}

template<class T> void verCatMatrix(vector< vector<T> > &res, vector < vector<T> > &A)
{
    if(res.size() == 0)
    {
       res = A;
    }
    else
    {
       res.insert(res.end(),A.begin(),A.end()); // concatenate
    }
}
