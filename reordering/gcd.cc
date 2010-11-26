#include <iostream>

using namespace std;

int main()
{
    int a(30);
    int b(50);

    int tmp;
    while (b)
    {
        tmp = b;
        b= a%b;
        a = tmp;
    }

    cout<<"gcd = "<<a<<endl;

    return 0;
}
