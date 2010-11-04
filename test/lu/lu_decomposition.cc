// Copyright (c) 2000-2008, Texas Engineering Experiment Station(TEES), a
// component of the Texas A&M University System.
// All rights reserved.
 
// Redistribution and use in source and binary forms, with or without
// modification, are not permitted without specific prior written permission
// from TEES.
 
// If written permission is obtained for redistribution or further use, the
// following conditions must be met:
 
// 1) Redistributions of source code must retain the above copyright notice,
// this list of conditions and the disclaimer below.
 
// 2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions, and the disclaimer below in the documentation and/or
// other materials provided with the distribution.
 
// 3) Neither the name of TEES, the name of the Texas A&M University System, nor
// the names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.
 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


#include "lu_decomposition.hh"
#include <iostream>

using namespace std;

LU::LU(int n)
{
    A.resize(n,vector<double> (n,1.0));
    pivot.resize(n,0);

    for (int i=0; i<n; i++)
        for (int j=n-2; j>=0; j--)
            A[i][j] = (i+1)*A[i][j+1];
}

void LU::decomposition()
{
    const int size(pivot.size());
    double max;

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

vector<double> LU::solve(vector<double> &b)
{
    // solve the linear equation Lx=b for x where L is a lower triangular
    // matrix
    const int size(pivot.size());
    vector<double> x(size,0.0);

    for (int k=0; k<size; k++)
    {
        if (pivot[k] !=k)
        {
            double tmp(b[k]);
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
            double tmp(b[k]);
            b[k]=b[pivot[k]];
            b[pivot[k]]=tmp;
        }
        for (int i=k+1; i<size; i++)
            x[k] -= x[i]*A[k][i];
    }

    return x;
}
