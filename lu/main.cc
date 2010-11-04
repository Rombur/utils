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
//

#include "crout_pivot.c"
#include "lu_decomposition.hh"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
//    double A[3][3];
//    double B[3];
//    double x[3];
//    int pivot[3];
//
//    A[0][0] = 1.0;
//    A[0][1] = 1.0;
//    A[0][2] = 1.0;
//    A[1][0] = 4.0;
//    A[1][1] = 2.0;
//    A[1][2] = 1.0;
//    A[2][0] = 9.0;
//    A[2][1] = 3.0;
//    A[2][2] = 1.0;
//
//    B[0] = 1.0;
//    B[1] = 0.0;
//    B[2] = 0.0;
//
//    for (int i=0; i<3;i++)
//    {
//        for (int j=0; j<3; j++)
//            cout<<A[i][j]<<"\t";
//        cout<<endl;
//    }         
//    int err;
//    err = Crout_LU_Decomposition_with_Pivoting(&A[0][0],pivot,3);
//    cout<<err<<endl;
//
//    for (int i=0; i<3;i++)
//    {
//        for (int j=0; j<3; j++)
//            cout<<A[i][j]<<"\t";
//        cout<<endl;
//    }
//    for (int i=0; i<3; i++)
//        cout<<pivot[i]<<"\t";
//    cout<<endl;
//
//    err = Crout_LU_with_Pivoting_Solve(&A[0][0],B,pivot,x,3);
//
//    for (int i=0; i<3; i++)
//        cout<<x[i]<<"\t";
//    cout<<endl;

    int n(3);
    LU lu(n);

    for (int i=0; i<n;i++)
    {
        for (int j=0; j<n; j++)
            cout<<lu.A[i][j]<<"\t";
        cout<<endl;
    }
    cout<<endl;
    lu.decomposition();
    cout<<setw(12);
    for (int i=0; i<n;i++)
    {
        for (int j=0; j<n; j++)
            cout<<lu.A[i][j]<<setw(12);
        cout<<endl;
    }
    cout<<endl;
    for (int i=0; i<n; i++)
        cout<<lu.pivot[i]<<"\t";
    cout<<endl;

    vector<double> x;
    vector<double> b;
    vector<vector<double> > C;
    for (int i=0; i<n; i++)
    {
        b.resize(n,0.0);
        b[i] = 1.0;
        x = lu.solve(b);
        C.push_back(x);
        b.clear();
    }

    cout<<endl;
    cout<<setw(12);
    for (int i=0; i<n;i++)
    {
        for (int j=0; j<n; j++)
            cout<<C[j][i]<<setw(12);
        cout<<endl;
    }
    cout<<endl;

    return 0;
}
