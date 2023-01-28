//bogo sort
#include<iostream>
#include<cstdlib>
#include<ctime>
using namespace std;

void printArray(int a[], int n) {
    for (int i=0; i<n; i++)
        cout<<" "<<a[i];
    cout<<endl;
}

bool isSorted(int a[], int n) {
    for (int i=0; i<n-1; i++)
        if (a[i]>a[i+1])
            return false;
    return true;
}

void shuffle(int a[], int n) { //Fisher Yates shuffling
    srand(time(0));
    for (int i=n-1; i>0; i--) {
        int j = rand()%(i+1); //j = a random integer with 0 <= j <= i
        swap(a[i], a[j]);
    }
}

void bogoSort(int a[], int n) {
    while (!isSorted(a, n))
        shuffle(a, n);
}

int main() {
    int a[] = {64, 34, 25, 12, 22, 11, 90};
    int n = sizeof(a)/sizeof(a[0]);
    bogoSort(a, n);
    cout<<"Sorted Array :";
    printArray(a, n);
    return 0;   
}