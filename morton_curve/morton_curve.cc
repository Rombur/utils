#include <iostream>

/* 
 * For 128 bits, there is not built-in type in C++ so this function need to be
 * commented.
 *
 * uint128_t separate_by_2(uint128_t u)
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

unsigned long long int separate_by_2(unsigned long long int u)
{
  // Keep the 21 least-significant bits
  u &= 0x1fffff;

  // Shift the bits and filter them
  u = (u | u<<32) & 0x1f00000000ffff;
  u = (u | u<<16) & 0x1f0000ff0000ff;
  u = (u | u<<8) & 0x100f00f00f00f00f;
  u = (u | u<<4) & 0x10c30c30c30c30c3;
  u = (u | u<<2) & 0x1249249249249249;

  return u;
}

unsigned long long int separate_by_1(unsigned long long int u)
{
  // Keep the 32 least-significant bits
  u &= 0xffffffff;

  // Shift the bits and filter them
  u = (u | u<<16) & 0x0000ffff0000ffff;
  u = (u | u<<8) & 0x00ff00ff00ff00ff;
  u = (u | u<<4) & 0x0f0f0f0f0f0f0f0f;
  u = (u | u<<2) & 0x3333333333333333;
  u = (u | u<<1) & 0x5555555555555555;

  return u;
}

// Return an unsigned long long int but only using 63 bits. The inputs are also unsigned
// long long int but we use only 21 bits.
unsigned long long int compute_morton_ordering(unsigned long long int x, 
    unsigned long long int y, unsigned long long int z)
{
  unsigned long long int splitted_x = separate_by_2(x);
  unsigned long long int splitted_y = separate_by_2(y);
  unsigned long long int splitted_z = separate_by_2(z);

  // Return the order
  return splitted_x | (splitted_y<<1) | (splitted_z<<2);
}

// Return an unsigned long long int. The inputs are also unsigned long long int
// but we only use 32 bits.
unsigned long long int compute_morton_ordering(unsigned long long int x,
    unsigned long long int y)
{
  unsigned long long int splitted_x = separate_by_1(x);
  unsigned long long int splitted_y = separate_by_1(y);

  // Return the order
  return splitted_x | (splitted_y<<1);
}

int main()
{
  unsigned long long int x1(5); 
  unsigned long long int y1(9); 
  unsigned long long int z1(1); 

  unsigned long long int x2(400);
  unsigned long long int y2(300);
  unsigned long long int z2(12);

  unsigned long long int x3(800);
  unsigned long long int y3(600);
  
  std::cout<<compute_morton_ordering(x1,y1,z1)<<" = "<<1095<<std::endl;
  std::cout<<compute_morton_ordering(x2,y2,z2)<<" = "<<52501888<<std::endl;

  std::cout<<compute_morton_ordering(x3,y3)<<" = "<<861824<<std::endl;

  return 0;
}
