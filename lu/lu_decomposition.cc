#include "lu_decomposition.hh"
#include <iostream>

using namespace std;

LU::LU(int n)
{
    A.resize(n,vector<long double> (n,1.0));
    pivot.resize(n,0);

    for (int i=0; i<n; i++)
        for (int j=n-2; j>=0; j--)
            A[i][j] = (i+1)*A[i][j+1];
}

void LU::decomposition()
{
    const int size(pivot.size());
    long double max;

    // for each row and column
    for (int k=0; k<size; k++)
    {
        // find the pivot row
        pivot[k]=k;
        max = fabs(A[k][k]);
        for (int j=k+1; j<size;j++)
            if (max < fabs(A[j][k]))
            {
                max = fabs(A[j][k]);
                pivot[k] = j;
            }
        // and if the pivot row differs from the current row, then interchange
        // the two rows
        if (pivot[k] != k)
        {
            const int piv(pivot[k]);
            for (int j=0; j<size; j++)
            {

               max = A[k][j];
               A[k][j] = A[piv][j];
               A[piv][j] = max;
            }
        }

        // find the upper triangular matrix elements for row k
        for (int j=k+1; j<size; j++)
            A[k][j] /= A[k][k];

        // update remaining matrix
        for (int i=k+1; i<size; i++)
            for (int j=k+1; j<size; j++)
                A[i][j] -= A[i][k]*A[k][j];
    }
}

vector<long double> LU::solve(vector<long double> &b)
{
    // solve the linear equation Lx=b for x where L is a lower triangular
    // matrix
    const int size(pivot.size());
    vector<long double> x(size,0.0);

    for (int k=0; k<size; k++)
    {
        if (pivot[k] !=k)
        {
            long double tmp(b[k]);
            b[k]=b[pivot[k]];
            b[pivot[k]]=tmp;
        }
        x[k] = b[k];
        for (int i=0; i<k; i++)
            x[k] -= x[i]*A[k][i];
        x[k] /=A[k][k];
    }
    // solve the linear equation Ux=y, where y is the solution obtained above
    // of Lx=b and U is an upper triangular matrix. The diagonal of the upper
    // triangular part of the matrix is assumed to be 1.
    for (int k=size-1; k >= 0; k--)
    {
        if (pivot[k] != k)
        {
            long double tmp(b[k]);
            b[k]=b[pivot[k]];
            b[pivot[k]]=tmp;
        }
        for (int i=k+1; i<size; i++)
            x[k] -= x[i]*A[k][i];
    }

    return x;
}
