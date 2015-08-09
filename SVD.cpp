#include<iostream>
#include<vector>
#include<limits>
#include<cmath>
#include<fstream>

//Created by tarmizi adam 9/8/2015
// This code perform Singular Valude Decomposition (SVD) and solves a linear system Ax = b
// using SVD. Where the matrix A is singular. This is done with the pesudoinverse.
// The code here is taken from the book:

//          Numerical Recipes: The art of scientific computing by Press, Teukolsky, Vetterling
//                  and Flannery. and: Numerical Recipes Webnote no 2., Rev 1 available at
//                  http://www.nr.com/webnotes/nr3web2.pdf

// From section 2.6 : Singular Value Decomposition. The codes here are just slightly modified from the original.


using namespace std;

template<class T>
inline T SIGN(const T &a, const T &b)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

template<class T>
inline const T &MAX(const T &a, const T &b)
{
    return b > a ? (b) : (a);
}

template<class T>
inline const T &MIN(const T &a, const T &b)
{
    return b < a ? (b) : (a);
}

template<class T>
inline T SQR(const T a)
{
    return a*a;
}

struct SVD
{
    int m;
    int n;

    vector< vector<double> > U; // Matrix U;
    vector< vector<double> > V; // Matrix V;

    vector<double> w;   // The diagonal matrix W

    double eps;
    double tsh;
    SVD(vector< vector<double> > &A);

    void solve(vector<double> &b, vector<double> &x, double thresh);
    void solve(vector< vector<double> > &b, vector< vector<double> > &x, double thresh);

    void decompose();
    void reorder();
    double pythag(const double a, const double b);

};

SVD::SVD(vector< vector<double> > &A) : m(A.size()),n(A[0].size()),U(A),V(n, vector<double>(n,0.0)),w(n)
{
    eps = numeric_limits<double>::epsilon();
    decompose();
    reorder();
    tsh = 0.5*sqrt(m+n+1.)*w[0]*eps;

}

void SVD::solve(vector<double> &b, vector<double> &x, double thresh)
{
    double s;

    if(b.size() != m || x.size() != n)
    {
        cout << "SVD solve bad sizes";
        return;
    }

    vector<double> tmp(n,0.0);

    tsh = (thresh >= 0. ? thresh: 0.5*sqrt(m+n+1.)*w[0]*eps);

    for(int j =0; j<n; j++)
    {
        s = 0.0;
        if(w[j] > tsh)
        {
            for(int i =0; i<m; i++)
            {
                s += U[i][j]*b[i];
            }
            s /= w[j];
        }

        tmp[j] =s;
    }

    for(int j =0; j<n; j++)
    {
        s = 0.0;

        for(int jj =0; jj<n; jj++)
        {
            s += V[j][jj]*tmp[jj];

        }

        x[j] = s;
    }
}

void SVD::solve(vector< vector<double> > &b, vector< vector<double> > &x, double thresh)
{
    int i;
    int j;
    int m = b[0].size(); // size of column

    if(b.size() != n || x.size() != n || b[0].size() != x[0].size())
    {
        cout << "SVD solve bad size !" << endl;
        return;
    }

    vector<double> xx(n,0.0);

    for(j =0; j<m; j++)
    {
        for(i =0; i<n; i++)
        {
            xx[i] = b[i][j];
            solve(xx,xx, thresh);

            for(i = 0; i<n; i++)
            {
                x[i][j] = xx[i];
            }
        }
    }
}

void SVD::decompose()
{
    bool flag;

    int i,its,j,jj,k,l,nm;

    double Anorm,c,f,g,h,s,scale,x,y,z;

    vector<double> rv1(n,0.0);


    for(i =0; i < n; i++)
    {
      l = i+2;
      rv1[i] = scale*g;
      g = 0.0;
      s = 0.0;
      scale = 0.0;

      if(i < m)
      {
          for(k = i; k < m; k++)
          {
              scale += fabs(U[k][i]);
          }
          if(scale != 0.0)
          {
              for(k = i; k<m; k++)
              {
                  U[k][i] /= scale;
                  s += U[k][i]*U[k][i];
              }

              f = U[i][i];
              g = -SIGN(sqrt(s),f);
              h = f*g-s;
              U[i][i] = f-g;

              for(j = l-1; j<n; j++)
              {
                  for(k = i, s = 0.0; k<m; k++)
                  {
                      s += U[k][i]*U[k][j];
                  }

                  f = s/h;

                  for(k = i; k<m; k++)
                  {
                     U[k][j] += f*U[k][i];
                  }

              }

              for(k = i; k<m; k++)
              {
                  U[k][i] *= scale;
              }
          }
      }

      w[i] = scale*g;
      g = 0.0;
      s = 0.0;
      scale = 0.0;

      if(i+1 <= m && i+1 != n)
      {
          for(k = l-1; k<n; k++)
          {
              scale += fabs(U[i][k]);
          }

          if(scale != 0.0)
          {
              for(k = l-1; k<n; k++)
              {
                  U[i][k] /= scale;
                  s += U[i][k]*U[i][k];
              }

              f = U[i][l-1];
              g = -SIGN(sqrt(s),f);
              h = f*g-s;
              U[i][l-1] = f-g;

              for(k = l-1; k<n; k++)
              {
                  rv1[k] = U[i][k]/h;
              }

              for(j = l-1; j<m; j++)
              {
                  for(s=0.0, k = l-1; k<n; k++)
                  {
                      s += U[j][k]*U[i][k];
                  }

                  for(k = l-1; k<n; k++)
                  {
                      U[j][k] += s*rv1[k];
                  }
              }

              for(k = l-1; k<n; k++)
              {
                  U[i][k] *= scale;
              }
          }
      }

        Anorm = MAX(Anorm,(fabs(w[i])+fabs(rv1[i])));

    }

    for(i = n-1; i>=0; i--)
    {
        if(i < n-1)
        {
            if(g != 0.0)
            {

                for(j = l; j<n; j++)
                {
                    V[j][i] = (U[i][j]/U[i][l])/g;
                }

                for(j =l; j<n; j++)
                {
                    for(s=0.0, k =l; k<n; k++)
                    {
                        s+=U[i][k]*V[k][j];
                    }

                    for(k =l; k<n; k++)
                    {
                        V[k][j] += s*V[k][i];
                    }
                }
            }

            for(j =l; j<n; j++)
            {
                V[i][j] = V[j][i] =0.0;
            }
        }

        V[i][i]=1.0;
        g=rv1[i];
        l=i;
    }

    for(i = MIN(m,n)- 1; i>=0; i--)
    {
        l = i+1;
        g = w[i];

        for(j =l; j<n; j++)
        {
            U[i][j] = 0.0;
        }

        if(g != 0.0)
        {
            g = 1.0/g;

            for(j = l; j<n; j++)
            {
                for(s=0.0, k =l; k<m; k++)
                {
                    s +=U[k][i]*U[k][j];
                }

                f = (s/U[i][i])*g;

                for(k =i; k<m; k++)
                {
                    U[k][j] += f*U[k][i];
                }
            }

            for(j =i; j<m; j++)
            {
                U[j][i] *= g;
            }
        }
        else
        {
            for(j  =i; j<m; j++)
            {
                U[j][i] =0.0;
            }
        }

        ++U[i][i];


    }

    for(k = n-1; k>=0; k--)
    {
        for(its = 0; its<30; its++)
        {
            flag = true;

            for(l = k; l>=0; l--)
            {
                nm = l-1;

                if(l == 0 || fabs(rv1[l]) <= eps*Anorm)
                {
                    flag = false;
                    break;
                }

                if(fabs(w[nm]) <= eps*Anorm)
                {
                    break;
                }
            }

            if(flag)
            {
                c = 0.0;
                s = 1.0;

                for(i = l; i<k+1; i++)
                {
                    f = s*rv1[i];
                    rv1[i] = c*rv1[i];

                    if(fabs(f) <= eps*Anorm)
                    {
                        break;
                    }

                    g = w[i];
                    h = pythag(f,g);
                    w[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;

                    for(j =0; j<m; j++)
                    {
                        y = U[j][nm];
                        z = U[j][i];
                        U[j][nm] = y*c+z*s;
                        U[j][i] = z*c-y*s;
                    }
                }
            }

            z = w[k];

            if(l == k)
            {
                if(z < 0.0)
                {
                    w[k] = -z;

                    for(j = 0; j<n; j++)
                    {
                        V[j][k] = -V[j][k];
                    }
                }

                break;
            }

            if(its == 29)
            {
                cout << "No convergence in 30 svdcmp iterations";
                return;
            }

            x = w[l];
            nm = k-1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f, 1.0);
            f = ((x-z)*(x+z)+ h*((y/(f+SIGN(g,f)))-h))/x;
            c = s = 1.0;

            for(j =l; j<=nm; j++)
            {
                i = j+1;
                g = rv1[i];
                y = w[i];
                h = s*g;
                g = c*g;
                z = pythag(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = g*c-x*s;
                h = y*s;
                y *= c;

                for(jj =0; jj<n; jj++)
                {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x*c+z*s;
                    V[jj][i] = z*c-x*s;
                }

                z = pythag(f,h);
                w[j] = z;

                if(z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }

                f = c*g+s*y;
                x = c*y-s*g;

                for(jj = 0; jj<m; jj++)
                {
                    y = U[jj][j];
                    z = U[jj][i];
                    U[jj][j] = y*c+z*s;
                    U[jj][i] = z*c-y*s;
                }
            }

            rv1[l] =0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

}

void SVD::reorder()
{
    int inc = 1;
    double sw;
    vector<double> su(m,0.0);
    vector<double> sv(n,0.0);

    do
    {
        inc *= 3;
        inc++;
    }while(inc <= n);

    do
    {
        inc /= 3;
        for(int i = inc; i<n; i++)
        {
            sw = w[i];

            for(int k  =0; k<m; k++)
            {
                su[k] = U[k][i];
            }
            for(int k =0; k<n; k++)
            {
                sv[k] = V[k][i];
            }

            int j = i;

            while(w[j-inc] < sw)
            {
                w[j] = w[j-inc];

                for(int k = 0; k<m; k++)
                {
                   U[k][j] = U[k][j-inc];
                }

                for(int k = 0; k<n; k++)
                {
                    V[k][j]=V[k][j-inc];
                }

                j -= inc;

                if(j < inc)
                {
                    break;
                }
            }

            w[j] = sw;

            for(int k = 0; k<m; k++)
            {
                U[k][j] = su[k];
            }

            for(int k =0; k<n; k++)
            {
                V[k][j] = sv[k];
            }
        }
    }while(inc > 1);

    for(int k =0; k<n; k++)
    {
        int s =0;

        for(int i =0; i<m; i++)
        {
            if(U[i][k] < 0.)
            {
                s++;
            }
        }

        for(int j =0; j<n; j++)
        {
            if(V[j][k]< 0.)
            {
                s++;
            }
        }

        if(s > (m+n)/2)
        {
            for(int i = 0; i<m; i++)
            {
                U[i][k] = -U[i][k];
            }

            for(int j = 0; j<n; j++)
            {
                V[j][k] = -V[j][k];
            }
        }
    }

    ofstream out("U.csv");

    for(int i =0; i < U.size(); i++)
    {
        for(int j =0; j<U[0].size();j++)
        {
            out << U[i][j] << ", ";
        }
        out << endl;

        //out << w[i] << ", ";
    }


}


double SVD::pythag(const double a, const double b)
{
    double absa = abs(a);
    double absb = abs(b);

    return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)): (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}

int main()
{
    ofstream sol("sol.txt");
    vector< vector<double> > a;
    vector<double> x(3,0.0);
    vector<double> b = {1,2,5,3};

    a = {{2.0, 2.0, 5.0},
         {4.0, 5.0, 1.0},
         {7.0, 8.0, 9.0},
         {13.0, 11.0, 12.0},
         };

    //a ={{2, 6},{1,3}}; // an example of singular matrix

    SVD s(a); // Decompose
    s.solve(b,x,-1.); // Solve the linear system Ax = b.

   for(int i = 0; i<x.size(); i++)
   {
        sol << x[i] << ", "; // Print out our solution x, to a txt file...
   }

    return 0;
}
