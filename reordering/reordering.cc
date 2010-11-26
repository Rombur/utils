#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;
typedef vector<d_vector> d_matrix;

void read(char* file, vector<d_matrix> &xs, d_vector &s, unsigned int &n_a, 
        unsigned int &n_b)
{
   unsigned int L_max;
   fstream input(file, ios::in);

   if (input)
   {
       input >> n_a >> n_b >> L_max;
       
       const unsigned int n(n_a+n_b);      
       s.resize(n,0.0);
       xs.resize(n);
       for (unsigned int i=0; i<n; i++)
           xs[i].resize(L_max,vector<double> (n,0.0));

       for (unsigned int i=0; i<n; i++)
           input >> s[i];

       for (unsigned int l=0; l<L_max; l++)
           for (unsigned int i=0; i<n; i++)
               for (unsigned int j=0; j<n; j++)
                   input >> xs[i][l][j];

       input.close();
   }
   else
       cerr << "Cannot open the file" << endl;
}

template <typename T>
void print_matrix(vector<T> const &xs)
{
    for (unsigned int j=0; j < xs[0].size(); j++)
    {
        for (unsigned int i=0; i < xs.size(); i++)
        {
            for (unsigned int k=0; k < xs[i][j].size(); k++)
                cout<<setw(3)<<xs[i][j][k];
            cout<<endl;
        }
        cout<<endl;
    }
}

unsigned int compute_gcd(unsigned int n_a, unsigned int n_b)
{
    unsigned int tmp;
    while (n_b)
    {
        tmp = n_b;
        n_b = n_a % n_b;
        n_a = tmp;
    }
    return n_a;
}

void create_linear_permutation(d_matrix &permutation_m, ui_vector &permutation_v,
        unsigned int const &n, unsigned int const &n_a, unsigned int const &n_b, 
        bool & identity)
{
    unsigned int gcd = compute_gcd(n_a,n_b);
    if (gcd != 1)
    {
        identity = false;
        unsigned int gcd_a = n_a/gcd;
        unsigned int gcd_b = n_b/gcd;

        for (unsigned int i=0; i<n_a; i++)
        {
            unsigned int column(0);

            column = i/gcd_a*(gcd_a+gcd_b)+i%gcd_a;
            permutation_m[i][column] = 1.0;
            permutation_v[i] = column;
        }

        for (unsigned int i=0; i<n_b; i++)
        {
            unsigned int column(gcd_a);

            column += i/gcd_b*(gcd_a+gcd_b)+i%gcd_b;
            permutation_m[i+n_a][column] = 1.0;
            permutation_v[i+n_a] = column;
        }
    }
    else
        for (unsigned int i=0; i<n; i++)
            permutation_v[i]=i;
    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
           cout<<setw(3)<<permutation_m[i][j];
        cout<<endl;
    }
    cout<<endl;
}

void matrix_multiplication(vector<d_matrix> &xs, d_matrix const &permutation_m, 
        const unsigned int &n)
{
    const unsigned int L_max(xs[0].size());
    vector<d_matrix> result;
    result.resize(n);
    for (unsigned int i=0; i<n; i++)
        result[i].resize(L_max, d_vector (n,0.0));
    
    // R = XS * P
    for (unsigned int i=0; i<n; i++)
        for (unsigned int j=0; j<L_max; j++)
            for (unsigned int k=0; k<n; k++)
                for (unsigned int m=0; m<n; m++)
                    result[i][j][m] += xs[i][j][k] * permutation_m[k][m];

    xs.clear();
    xs.resize(n);
    for (unsigned int i=0; i<n; i++)
        xs[i].resize(L_max, vector<double> (n,0.0));

    // R = P^t * XS
    for (unsigned int i=0; i<n; i++)
        for (unsigned int m=0; m<n; m++)
            for (unsigned int j=0; j<L_max; j++)
                for (unsigned int k=0; k<n; k++)
                    xs[i][j][k] += permutation_m[m][i] * result[m][j][k];
}

d_matrix reordering_linear(vector<d_matrix> &xs, ui_vector &permutation_v,
        unsigned int const &n_a,unsigned int const &n_b)
{
    bool identity(true);
    const unsigned int n(n_a+n_b);
    d_matrix permutation_m(n,vector<double> (n,0.0));
    permutation_v.resize(n);

    create_linear_permutation(permutation_m,permutation_v,n,n_a,n_b,identity);
    
    if (!identity)
        matrix_multiplication(xs,permutation_m,n);

    return permutation_m;
}

void create_log_permutation(d_matrix &permutation_m, ui_vector &permutation_v,
        unsigned int const &n, unsigned int const &n_a, unsigned int const &n_b)
{
    for (unsigned int i=0; i<n_a; i++)
    {
        if (i<n_b)
        {
            permutation_m[i][2*i] = 1.0;
            permutation_v[i] = 2*i;
        }
        else
        {
            permutation_m[i][2*n_b+i] = 1.0;
            permutation_v[i] = 2*n_b+i;
        }
    }

    for (unsigned int i=0; i<n_b; i++)
    {
        if (i<n_a)
        {
            permutation_m[i+n_a][2*i+1] = 1.0;
            permutation_v[i+n_a] = 2*i+1;
        }
        else
        {
            permutation_m[i+n_a][n_a+i] = 1.0;
            permutation_v[i+n_a] = n_a+i;
        }
    }

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
           cout<<setw(3)<<permutation_m[i][j];
        cout<<endl;
    }
    cout<<endl;
}

d_matrix reordering_log(vector<d_matrix> &xs, ui_vector &permutation_v,
        unsigned int const &n_a, unsigned int const &n_b)
{
    const unsigned int n(n_a+n_b);
    d_matrix permutation_m(n,d_vector (n,0.0));
    permutation_v.resize(n);

    create_log_permutation(permutation_m,permutation_v,n,n_a,n_b);
    
    matrix_multiplication(xs,permutation_m,n);

    return permutation_m;
}

template <typename T>
void print_vector(T const &s)
{
    const unsigned int i_max(s.size());
    for (unsigned int i=0; i<i_max; i++)
        cout<<s[i]<<endl;
    cout<<endl;
}

void mv_multiplication(ui_vector const &permutation_v, d_vector &s)
{
    d_vector tmp(s);
    const unsigned int n(s.size());
    
    for (unsigned int i=0; i<n; i++)
        s[permutation_v[i]] = tmp[i];
}

int main(int argc, char **argv)
{
    unsigned int n_energy_groups_a(0), n_energy_energy_groups_b(0);
    ui_vector permutation_v;
    d_vector s;
    d_matrix permutation;
    vector<d_matrix> xs;
    
    read(argv[2],xs,s,n_energy_groups_a,n_energy_energy_groups_b);

    print_matrix(xs);

    if (strcmp(argv[1],(char*)"linear") == 0)
        permutation = reordering_linear(xs,permutation_v,n_energy_groups_a, 
                n_energy_energy_groups_b);
    else
        permutation = reordering_log(xs,permutation_v,n_energy_groups_a, 
                n_energy_energy_groups_b);

    print_matrix(xs);

    cout<<"--------------------------------------\n";

    print_vector(s);

    mv_multiplication(permutation_v,s);

    print_vector(permutation_v);

    print_vector(s);

    return 0;
}
