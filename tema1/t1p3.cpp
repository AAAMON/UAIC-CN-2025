// Task: Polinome hierarchy based on approximations of the sine function

#include <iostream>

// constants:
double c1 = 0.16666666666666666666666666666667;
double c2 = 0.00833333333333333333333333333333;
double c3 = 1.984126984126984126984126984127e-4;
double c4 = 2.7557319223985890652557319223986e-6;
double c5 = 2.5052108385441718775052108385442e-8;
double c6 = 1.6059043836821614599392377170155e-10;

double P1(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + c2 * y));
}

double P2(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 - c3 * y)));
}

double P3(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 + y * (- c3  + c4 * y))));
}

double P4(double x)
{
    double y = x * x;
    return x * (1 + y * (-0.166 + y * (0.00833 + y * (- c3  + c4 * y))));
}

double P5(double x)
{
    double y = x * x;
    return x * (1 + y * (-0.1666 + y * (0.008333 + y * (- c3  + c4 * y))));
}

double P6(double x)
{
    double y = x * x;
    return x * (1 + y * (-0.16666 + y * (0.0083333 + y * (- c3  + c4 * y))));
}

double P7(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 + y * (- c3  + y * (c4 - c5 * y)))));
}

double P8(double x)
{
    double y = x * x;
    return x * (1 + y * (- c1 + y * (c2 + y * (- c3  + y * (c4 + y * (- c5 + c6 * y))))));
}

int main ()
{
    // something

    return 0;
}
