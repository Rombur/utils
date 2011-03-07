#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

enum type 
{
  centers, lines, nodes
};

// number as to be a float, double or long double
template <typename number>
class lipschitz
{
  public :
    lipschitz(int lvl, int dim);
    void mapd(number const &x, vector<number> &y, type key);
    void node(number i);
    void invmad(vector<number> & xp, int & kxx, vector<number> &p, 
        int incr);
    void xyd(number &xx, vector<number> &y);
    void numbr(int &iss);

  private :
    int level;
    int dimension;
    int l;
    int n_exp;
    int iq;   // unsigned int ?
    number del;
    vector<number> iu;
    vector<number> iv;
};

template <typename number>
lipschitz<number>::lipschitz(int lvl, int dim) : 
  level(lvl),
  dimension(dim),
  iu(dim,0.0),
  iv(dim,0.0)
  {}

template <typename number>
void lipschitz<number>::mapd(number const &x, vector<number> &y, type key)
{
  int n1 = dimension-1;
  number p(0.0);
  //n_exp = 2^dimension
  n_exp = pow(2.0,dimension);
  number d(x);
  number r=0.5;
  int it(0);
  number dr(n_exp);
  // mne = dr^level
  number mne(pow(dr,level));
  vector<int> iw(dimension,1);
  y.clear();
  y.resize(dimension,0.0);
  int k(0);
  number dd;
  int is;
  int i; // int PAS unsigned

  switch (key) 
  {
    case centers :
      break;

    case lines :
      d = d*(1.0-1.0/mne);
      k = 0;
      break;

    case nodes :
      dr = mne/n_exp;
      dr = dr - fmod(dr,1.0);
      dd = mne - dr;
      dr = d*dd;
      dd = dr - fmod(dr,1.0);
      dr = dd + (dd-1.0)/(n_exp-1.0);
      dd = dr - fmod(dr,1.0);
      d = dd*(1./mne);
  }

  for (int j=0; j<level; ++j)
  {
    iq = 0.0;
    if (x == 1.0)
    {
      is = n_exp - 1;
      d = 0.0;
    }
    else
    {
      d = d*n_exp;
      is = d;
      d = d - is;
    }

    i = is;
    node(i);

    i=iu[0];
    iu[0] = iu[it];
    iu[it] = i;

    i=iv[0];
    iv[0] = iv[it];
    iv[it] = i;

    if (l == 0)
      l = it;
    else if (l==it)
      l=0;

    if ( (iq>0)||((iq==0)&&(is==0)))
      k=1;
    else if (iq<0)
      k = (it == n1) ? 0 : n1;

    r = r/2.0;
    it = l;
    for (int kk=0; kk<dimension; ++kk)
    {
      iu[kk] = iu[kk]*iw[kk];
      iw[kk] = -iv[kk]*iw[kk];
      p = r*iu[kk];
      p = p + y[kk];
      y[kk] = p;
    }
  }

  if (key == lines)
  {
    if ( is == (n_exp-1) )
      i = -1;
    else 
      i = 1;
    p = 2.0*i*iu[k]*r*d;
    p = y[k] - p;
    y[k] = p;
  }
  else if (key == nodes)
  {
    for (int kk=0; kk<dimension; ++kk)
    {
      p = r*iu[kk];
      p += y[kk];
      y[kk] = p;
    }
  }
}

template<typename number>
void lipschitz<number>::node(number is)
{
  int iff;
  int k1,k2,j;

  if (is == 0)
  {
    l = dimension - 1.0;
    for (int i=0; i<dimension; ++i)
    {
      iu[i] = -1.0;
      iv[i] = -1.0;
    }
  }
  else if (is == n_exp-1.0)
  {
    l = dimension - 1.0;
    iu[0] = 1.0;
    iv[0] = 1.0;
    for (int i=1; i<dimension; ++i)
    {
      iu[i] = -1.0;
      iv[i] = -1.0;
    }
    iv[dimension-1] = 1.0;
  }
  else
  {
    iff = n_exp;
    k1 = -1;
    for (int i=0; i<dimension; ++i)
    {
      iff /= 2.0;
      if (is >= iff)
      {
         if ( (is==iff) && (is != 1))
         {
           l =i;
           iq =-1.0;
         }
         is -= iff;
         k2 =1;
      }
      else
      {
        k2 =-1;
        if ( (is==(iff-1.0)) && (is !=0) )
        {
          l = i;
          iq = 1.0;
        }
      }
      j=-k1*k2;
      iv[i] = j;
      iu[i] = j;
      k1 = k2;
    }
    iv[l] = iv[l]*iq;
    iv[dimension-1] = -iv[dimension-1];
  }
}

template <typename number>
void lipschitz<number>::invmad(vector<number> & xp, int & kxx, 
    vector<number> &p, int incr)
{
  /**
   * Preimages calculation :
   *  - xp = preimages ti be calculated
   *  - kxx = number of preimages being caculated
   *  - p = images for which preimages are calculated
   *  - incr = minimum number of map nodes that must bebetween preimages
   */ 
  number mne, d1, dd, x, dr;
  double r,d;
  vector<double> u(dimension,-1.0);
  vector<double> y(dimension,0.0);
  int kx;
  int k;
  int i;
  int kp = xp.size()-1;
  kx = 0;
  n_exp = pow(2,dimension);
  dr = n_exp;
  mne =1;
  r =0.5;
  for (i=0; i<level; ++i)
  {
    mne *= dr;
    r *= 0.5;
  }
  dr = mne/n_exp;
  dr = dr -fmod(dr,1.0);
  del = 1./(mne-dr);
  d1 = del*(incr+0.5);

  kx =-1;
  while (kx < kp)
  {
    for (i=0; i<dimension; ++i)
    {
      d = p[i];
      y[i] = d -r*u[i];
    }
    for (i=0; (i<dimension) && (fabs(y[i]) < 0.5); ++i); // ????
    if (i>= dimension) 
    {
      xyd(x,y);
      dr =x*mne;
      dd = dr- fmod(dr,1.0);
      dr = dd /n_exp;
      dd = dd -dr + fmod(dr,1.0);
      x = dd*del;
      if (kx > kp) 
        break;
      k=kx++;
      if (kx == 0)
        xp[0] = x;
      else
      {
        while (k>=0)
        {
          dr = fabs(x-xp[k]);
          if (dr <= d1)
          {
            for (kx--; k<kx; k++,xp[k]=xp[k+1]);
            goto m6;
          }
          else
            if (x<=xp[k])
            {
              xp[k+1] = xp[k];
              k--;
            }
            else
              break;
        }
        xp[k+1] = x;
      }
    }
m6 : for (i=dimension-1; (i>=0) && (u[i]=(u[i]<=0.0) ? 1 : -1)<0; --i);
     if (i<0)
       break;
  }
  kxx = ++kx;
}

template <typename number>
void lipschitz<number>::xyd(number &xx, vector<number> &y)
{
  /**
   * Calcuate preimage xx for the nearest m-level center of y (xx - left boundary 
   * point of m-level interval)
   */

  number x, r1;
  double r;
  vector<int> iw(10,1);
  int i,j,it,is;

  r =0.5;
  r1 =1.0;
  x = 0.0;
  it =0;
  for (j=0; j<level; j++)
  {
    r*=0.5;
    for (i=0; i<dimension; ++i)
    {
      iu[i] = (y[i]<0) ? -1 : 1;
      y[i] -=r*iu[i];
      iu[i] *=iw[i];
    }
    i = iu[0];
    iu[0] = iu[it];
    iu[it] =i;
    numbr(is);
    i=iv[0];
    iv[0] = iv[it];
    iv[it] = i;
    for (i=0; i<dimension; ++i)
      iw[i] = -iw[i]*iv[i];
    if (l==0)
      l=it;
    else if (l==it)
      l=0;
    it =l;
    r1=r1/n_exp;
    x+=r1*is;
  }
  xx=x;
}

template <typename number>
void lipschitz<number>::numbr(int &iss)
{
  /**
   * calculate s(u) = is, l(u), v(u) =iv by u=iu
   */
  int i,is,iff,k1,k2,l1;

  is =0;
  k1=-1;
  iff = n_exp;
  for (i=0; i<dimension; ++i)
  {
    iff = iff/2;
    k2 =-k1*iu[i];
    iv[i] = iu[i];
    k1=k2;
    if (k2<0)
      l1=i;
    else
    {
      is += iff;
      l=i;
    }
  }
  if (is==0)
    l=dimension-1;
  else
  {
    iv[dimension-1] = -iv[dimension-1];
    if (is==(n_exp-1))
      l=dimension-1;
    else if (l1==dimension-1)
      iv[l]=-iv[l];
    else
      l=l1;
  }
  iss = is;
}


int main(int argc, char* argv[])
{

  int n=2;
  int m=10;
  vector<double> y;
  double x(0.5);
  lipschitz<double> l(m,n);
  type t(nodes);

  l.mapd(x,y,t);

  for (int i=0; i<n; ++i)
    cout<<y[i]<<endl;
  cout<<endl;

  vector<double> xp(4,0.0);

  int kx = 4;
  l.invmad(xp,kx,y,1);

  cout<<"kx "<<kx<<endl;
  vector<double>::iterator end(xp.end());
  vector<double>::iterator xp_val(xp.begin());
  for ( ; xp_val<end; xp_val++)
  {
    cout<<*xp_val<<endl;
  }

  cout<<"------------------------------"<<endl;
  m = 14;
  n = 3;
  x = 0.55;

  lipschitz<double> l2(m,n);
  l2.mapd(x,y,t);

  for (int i=0; i<n; ++i)
    cout<<y[i]<<endl;

  xp.clear();
  xp.resize(8,0.0);
  kx = 8;
  l2.invmad(xp,kx,y,1);

  cout<<"kx "<<kx<<endl;
  end = xp.end();
  xp_val = xp.begin();
  for ( ; xp_val<end; xp_val++)
  {
    cout<<*xp_val<<endl;
  }

  return 0;
}                                                                            
