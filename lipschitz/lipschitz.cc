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
    lipschitz(unsigned int lvl, unsigned int dim, type KEY);
    void mapd(number const &x, vector<number> &y);
    void invmad(vector<number> & xp, int & kxx, vector<number> &p, 
        int incr);

  private :
    const unsigned int level;
    const unsigned int dimension;
    unsigned int l;
    unsigned int n_exp;
    int iq;   
    type key;
    number del;
    vector<number> iu;
    vector<number> iv;

    void node(number i);
    void xyd(number &xx, vector<number> &y);
    void numbr(int &iss);
};

template <typename number>
lipschitz<number>::lipschitz(unsigned int lvl, unsigned int dim, type KEY) : 
  level(lvl),
  dimension(dim),
  n_exp(pow(2,dim)),
  key(KEY),
  iu(dim,0.0),
  iv(dim,0.0)
  {}

template <typename number>
void lipschitz<number>::mapd(number const &x, vector<number> &y)
{
  unsigned int n1(dimension-1);
  unsigned int it(0);
  unsigned int k(0);
  unsigned int is;
  int i; 
  number p(0.0);
  number d(x);
  number r=0.5;
  number dr(n_exp);
  number mne(pow(dr,level));
  number dd;
  const number one(1.0);
  vector<int> iw(dimension,1);

  y.clear();
  y.resize(dimension,0.0);

  switch (key) 
  {
    case centers :
      break;

    case lines :
      d = d*(one-one/mne);
      k = 0;
      break;

    case nodes :
      dr = mne/n_exp;
      dr = dr - fmod(dr,one);
      dd = mne - dr;
      dr = d*dd;
      dd = dr - fmod(dr,one);
      dr = dd + (dd-one)/(n_exp-one);
      dd = dr - fmod(dr,one);
      d = dd*(one/mne);
  }

  for (unsigned int j=0; j<level; ++j)
  {
    iq = 0;
    if (x == one)
    {
      is = n_exp - 1;
      d = one;
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
    for (unsigned int kk=0; kk<dimension; ++kk)
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
    for (unsigned int kk=0; kk<dimension; ++kk)
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
  int iff,k1,k2,j;

  if (is == 0)
  {
    l = dimension - 1.0;
    for (unsigned int i=0; i<dimension; ++i)
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
    for (unsigned int i=1; i<dimension; ++i)
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
    for (unsigned int i=0; i<dimension; ++i)
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
  int kx(0);
  int k;
  int i;
  int dim(dimension);
  int kp = xp.size()-1;
  const number one(1.0);
  number d1,dd,x,d;
  number dr(n_exp);
  number mne(one);
  number r(0.5);
  vector<number> u(dimension,-1.0);
  vector<number> y(dimension,0.0);

  for (i=0; i<int(level); ++i)
  {
    mne *= dr;
    r *= 0.5;
  }
  dr = mne/n_exp;
  dr = dr -fmod(dr,one);
  del = one/(mne-dr);
  d1 = del*(incr+0.5);

  kx =-1;
  while (kx < kp)
  {
    for (unsigned int ii=0; ii<dimension; ++ii)
    {
      d = p[ii];
      y[ii] = d -r*u[ii];
    }
    for (i=0; (i<dim) && (fabs(y[i]) < 0.5); ++i);
    if (i>= dim) 
    {
      xyd(x,y);
      dr =x*mne;
      dd = dr- fmod(dr,one);
      dr = dd /n_exp;
      dd = dd -dr + fmod(dr,one);
      x = dd*del;
      if (kx > kp) 
        break;
      k=kx++;
      if (kx == 0)
        xp[0] = x;
      else
      {
        bool m6(false);
        while (k>=0 && !m6)
        {
          dr = fabs(x-xp[k]);
          if (dr <= d1)
          {
            for (kx--; k<kx; k++)
              xp[k+1] = xp[k+2];
            m6 = true;
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
        if(!m6)
          xp[k+1] = x;
      }
    }
    for (i=dimension-1; (i>=0) && (u[i]=(u[i]<=0.0) ? 1 : -1)<0; --i);
    if (i<0)
      break;
  }
  kxx = ++kx;
}

template <typename number>
void lipschitz<number>::xyd(number &xx, vector<number> &y)
{
  /**
   * Calculate preimage xx for the nearest m-level center of y (xx - left boundary 
   * point of m-level interval)
   */

  unsigned int it(0);
  int i,is;
  number r(0.5);
  number r1(1.0);
  number x(0.0);
  vector<int> iw(dimension,1);

  for (unsigned int j=0; j<level; j++)
  {
    r*=0.5;
    for (unsigned int ii=0; ii<dimension; ++ii)
    {
      iu[ii] = (y[ii]<0) ? -1 : 1;
      y[ii] -=r*iu[ii];
      iu[ii] *=iw[ii];
    }
    i = iu[0];
    iu[0] = iu[it];
    iu[it] =i;
    numbr(is);
    i=iv[0];
    iv[0] = iv[it];
    iv[it] = i;
    for (unsigned int ii=0; ii<dimension; ++ii)
      iw[ii] = -iw[ii]*iv[ii];
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
  unsigned int l1;
  unsigned int is(0);
  int k2;
  int k1(-1);
  int iff(n_exp);

  for (unsigned int i=0; i<dimension; ++i)
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

  unsigned int n=2;
  unsigned int m=10;
  vector<long double> y;
  long double x(0.5);
  type t(nodes);
  lipschitz<long double> l(m,n,t);

  l.mapd(x,y);

  for (unsigned int i=0; i<n; ++i)
    cout<<y[i]<<endl;
  cout<<endl;

  vector<long double> xp(4,0.0);

  int kx = 4;
  l.invmad(xp,kx,y,1);

  cout<<"kx "<<kx<<endl;
  vector<long double>::iterator end(xp.end());
  vector<long double>::iterator xp_val(xp.begin());
  for ( ; xp_val<end; xp_val++)
  {
    cout<<*xp_val<<endl;
  }

  cout<<"------------------------------"<<endl;
  m = 14;
  n = 3;
  x = 0.55;

  lipschitz<long double> l2(m,n,t);
  l2.mapd(x,y);

  for (unsigned int i=0; i<n; ++i)
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
