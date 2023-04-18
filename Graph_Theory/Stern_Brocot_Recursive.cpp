#include <cstdio>
#include <string>
#include <iostream>

using namespace std;

int a, b;
string str;
void Stern_Brocot(int num1, int den1, int num2, int den2) {
    if (a == num1 + num2 && b == den1 + den2)
        cout << str << endl;
    else if (a * (den1 + den2) > b * (num1 + num2)) {
        str += 'R';
        Stern_Brocot(num1 + num2, den1 + den2, num2, den2);
    } else {
        str += 'L';
        Stern_Brocot(num1, den1, num1 + num2, den1 + den2);
    }
}

int main() {
    scanf("%d %d", &a, &b);
    if (a == 1 && b == 1)
        return 0;

    int num1 = 0, den1 = 1, num2 = 1, den2 = 0;
    str = "";
    Stern_Brocot(num1, den1, num2, den2);
}