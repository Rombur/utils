#include <cmath>
#include <vector>

using namespace std;

class LU
{
    public :
        
        LU(int n);

        void decomposition();

        vector<long double> solve(vector<long double> &b);

        vector<vector<long double> > A;
        vector<int> pivot;

};

