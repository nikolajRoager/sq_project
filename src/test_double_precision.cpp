#include<iostream>
#include<cmath>
#include<float.h>

using namespace std;

int main()
{

    double x=1;
    double px=x;

    //Just divide by 2 until adding this to 1 does nothing
    while(1.0+x!=1.0)
    {
        px=x;//Previous x, which did still register when added to 1
        x/=2.0;
    }

    cout<<"Precision limit of double "<<px<<endl;
    cout<<"Official limit is 2^"<<DBL_MANT_DIG<<endl;
}
