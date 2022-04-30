#include <stdio.h>
#include <algorithm>
#include <vector>
#include <math.h>
#define COM 0
#define PRI 1
using namespace std;
vector<bool> a(110000001, 0);

int primes[9000001];
void sieve()
{
    int i, j, ctr;
    for (i = 6; i <= 100000000; i += 6) {
        a[i + 1] = PRI;
        a[i + 5] = PRI;
    }

    a[0] = a[1] = COM;
    a[2] = a[3] = a[5] = PRI;
    a[25] = COM;
    for (j = 35; j < 100000000; j += 30)
        a[j] = a[j + 20] = COM;
    int t;
    for (i = 6; i < 10000; i += 6) {
        if (a[i + 1] == PRI)
            for (j = (i + 1) * (i + 1); j < 100000000; j += 6 * (i + 1)) {
                a[j] = a[j + 4 * (i + 1)] = COM; //printf("%d ",j);scanf("%d",&t);
            }
        if (a[i + 5] == PRI)
            for (j = (i + 5) * (i + 5); j < 100000000; j += 6 * (i + 5)) {
                a[j] = a[j + 2 * (i + 5)] = COM; //printf("%d ",j);scanf("%d",&t);
            }
    }

    primes[0] = 2;
    primes[1] = 3;
    primes[2] = 5;
    ctr = 2;

    for (i = 6; i < 100000000; i += 6) {
        if (a[i + 1] == PRI) {
            ctr++;
            primes[ctr] = i + 1;
        }

        if (a[i + 5] == PRI) {
            ctr++;
            primes[ctr] = i + 5;
        }
    }
    return;
}

int main()
{
    segmented_sieve();
    for (int i = 0; i < 100; i += 100)
        printf("%d\n", primes[i]);
    return 0;
}
