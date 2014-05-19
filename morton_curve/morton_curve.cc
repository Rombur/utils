#include <iostream>

/* 
 * For 128 bits, there is not built-in type in C++ so this function need to be
 * commented.
 *
 * uint128_t split_by_3bits(uint128_t u)
 * {
 *   // Keep the 42 least-signigicant bits
 *   u &= 0x3ffffffffff;
 *
 *   u = (u | u<<64) & 0x3ff0000000000000000ffffffff;
 *   u = (u | u<<32) & 0x3ff00000000ffff00000000ffff;
 *   u = (u | u<<16) & 0x30000ff0000ff0000ff0000ff0000ff;
 *   u = (u | u<<8) & 0x300f00f00f00f00f00f00f00f00f00f;
 *   u = (u | u<<4) & 0x30c30c30c30c30c30c30c30c30c30c3;
 *   u = (u | u<<2) & 0x9249249249249249249249249249249;
 *
 *   return u;
 * }
 */

unsigned long long int split_by_3bits(unsigned long long int u)
{
  // Keep the 21 least-significant bits
  u &= 0x1fffff;

  u = (u | u<<32) & 0x1f00000000ffff;
  u = (u | u<<16) & 0x1f0000ff0000ff;
  u = (u | u<<8) & 0x100f00f00f00f00f;
  u = (u | u<<4) & 0x10c30c30c30c30c3;
  u = (u | u<<2) & 0x1249249249249249;

  return u;
}


// Return an unsigned long long int but only using 63 bits. The inputs are also unsigned
// long long int but we use only 21 bits.
unsigned long long int compute_morton_ordering(unsigned long long int x, 
    unsigned long long int y, unsigned long long int z)
{
  unsigned long long int splitted_x = split_by_3bits(x);
  unsigned long long int splitted_y = split_by_3bits(y);
  unsigned long long int splitted_z = split_by_3bits(z);

  // Return the order
  return splitted_x | (splitted_y<<1) | (splitted_z<<2);
}

int main()
{
  unsigned long long x(5); 
  unsigned long long y(9); 
  unsigned long long z(1); 
  
  std::cout<<compute_morton_ordering(x,y,z)<<std::endl;

  return 0;
}
