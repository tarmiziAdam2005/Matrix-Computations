#include <iostream>
#include <vector>

//Created 14/8/2015 by Tarmizi Adam. Simple program to create an Identity matrix

using namespace std;

vector< vector<double> > eye(size_t N);

int main()
{
    size_t N = 8;

    vector< vector<double> > I = eye(N);



    for(size_t row = 0; row < I.size(); row++)
    {
        for(size_t col =0; col < I[0].size(); col++)
        {
            cout << I[row][col] << ",";
        }
        cout << endl;
    }


    return 0;
}

vector< vector<double> > eye(size_t N)
{
    vector < vector<double> > identity(N,vector<double>(N,0.0));

    for(size_t row =0; row < N ; row++)
    {
        identity[row][row] = 1.0;
    }

    return identity;

}
