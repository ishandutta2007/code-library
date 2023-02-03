/*
Prime counting function implementation based on Linnik's Identity, and the
Dirichlet Hyperbola Method.

Nathan McKenzie, 5-18-2012
nathan _AT_ icecreambreakfast.com

*/

#include <bits/stdc++.h>
using namespace std;

typedef long long BigInt;
const int g_wheelLargestPrime = 19;
int scaleNum = 10;

void MakeWheel(int);
BigInt primeCountingFunction(BigInt);

int main(int argc, char *argv[]) {
  MakeWheel(g_wheelLargestPrime);

  int oldClock = (int)clock();
  int lastDif = 0;

  printf(
      "                                                                    ");
  printf("Time\n");
  printf(
      "                                                                    ");
  printf("Increase\n");
  printf(
      "                                                                    ");
  printf("for x%d\n", scaleNum);
  printf("         ");
  printf("__ Input Number __   __ Output Number __ _ MSec _ _ Sec _  Input\n");
  printf(
      "                                                                    \n");
  for (BigInt i = scaleNum; i <= 1000000000000000000; i *= scaleNum) {
    // printf( "%17I64d(10^%4.1f): ", i, log( (double)i )/log(10.0) );
    BigInt total = (BigInt)(primeCountingFunction(i) + 0.00001);
    int newClock = (int)clock();
    // printf( " %20I64d %8d : %4d: x%f\n",
    //     total, newClock - oldClock, ( newClock - oldClock ) / CLK_TCK,
    //     ( lastDif ) ? (double)( newClock - oldClock ) / (double)lastDif : 0.0
    //     );
    lastDif = newClock - oldClock;
    oldClock = newClock;
  }

  // getch();

  return 0;
}

/* ________________________________________________________________________________________

This is a pretty fast prime counting algorithm considering it uses O(episilon)
memory.
I've tried, from a number of different angles, to come up with ways of trading
off
memory usage to increase its runtime performance, but I haven't found any such
approaches.

I'm not entirely sure what the Big O runtime performance of this algorithm is,
ultimately... Up to around 10^15, it seems to be running in something like
O(n^4/5),
but it also seems to be taking slightly longer each power of 10.
________________________________________________________________________________________
*/

const double EPSILON = .00000000001;

typedef long long BigInt;

/* A wheel including 19 has 9.6 million entries. */

/* Used for the construction of the wheel - include more primes as needed,
but this is already enough primes to consume over 75 gigs of RAM */

const int g_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};

int g_wheelCycleEntries;
int g_wheelCyclePeriod;
int g_wheelFirstPrime;
int g_wheelBasePrimes;
int *g_wheelTranslation = 0;
int *g_wheelOffsets = 0;

BigInt g_latticepoints;
BigInt g_minVarValue;
BigInt g_boundary;

BigInt g_scale;
BigInt g_divisor;

BigInt g_lastMax;

int g_variablesLeft;
BigInt g_lastScaleDivisor;
BigInt g_scaleVal;

int g_mu[] = {
    0, 1,  -1, -1, 0, -1, 1,  -1, 0, 0,  1,  -1, 0, -1, 1,  1,  0, -1, 0,  -1,
    0, 1,  1,  -1, 0, 0,  1,  0,  0, -1, -1, -1, 0, 1,  1,  1,  0, -1, 1,  1,
    0, -1, -1, -1, 0, 0,  1,  -1, 0, 0,  0,  1,  0, -1, 0,  1,  0, 1,  1,  -1,
    0, -1, 1,  0,  0, 1,  -1, -1, 0, 1,  -1, -1, 0, -1, 1,  0,  0, 1,  -1, -1,
    0, 0,  1,  -1, 0, 1,  1,  1,  0, -1, 0,  1,  0, 1,  1,  1,  0, -1, 0,  0,
    0, -1, -1, -1, 0, -1, 1,  -1, 0, -1, -1, 1,  0, -1, -1, 1,  0, 0,  1,  1,
    0, 0,  1,  1,  0, 0,  0,  -1, 0, 1,  -1, -1, 0, 1,  1,  0,  0, -1, -1, -1,
    0, 1,  1,  1,  0, 1,  1,  0,  0, -1, 0,  -1, 0, 0,  -1, 1,  0, -1, 1,  1,
    0, 1,  0,  -1, 0, -1, 1,  -1, 0, 0,  -1, 0,  0, -1, -1, 0,  0, 1,  1,  -1,
    0, -1, -1, 1,  0, 1,  -1, 1,  0, 0,  -1, -1, 0, -1, 1,  -1, 0, -1, 0,  -1,
    0, 1,  1,  1,  0, 1,  1,  0,  0, 1,  1,  -1, 0, 1,  1,  1,  0, 1,  1,  1,
    0, 1,  -1, -1, 0, 0,  1,  -1, 0, -1, -1, -1, 0, -1, 0,  1,  0, 1,  -1, -1,
    0, -1, 0,  0,  0, 0,  -1, 1,  0, 1,  0,  -1, 0, 1,  1,  -1, 0,
};

/* Note that with 64 bit ints, we can't go above factorial( 20 ) anyway. */
BigInt g_factorial[] = {0,
                        1,
                        2,
                        6,
                        24,
                        120,
                        720,
                        5040,
                        40320,
                        362880,
                        3628800,
                        39916800,
                        479001600,
                        6227020800,
                        87178291200,
                        1307674368000,
                        20922789888000,
                        355687428096000,
                        6402373705728000,
                        121645100408832000,
                        2432902008176640000};

inline BigInt InversePower(BigInt x, BigInt y) {
  return ((BigInt)(pow((double)x + EPSILON, (1.0 / (double)y)) + EPSILON));
}

inline int GetwheelCyclePeriod(int cap) {
  int val = 1;
  int i = 0;

  while (g_primes[i] <= cap) {
    val *= g_primes[i];
    i++;
  }
  return val;
}

inline int GetFirstIncludedPrime(int cap) {
  int i = 0;

  while (g_primes[i] <= cap) {
    i++;
  }
  return g_primes[i];
}

inline int GetBasePrimes(int cap) {
  int i = 0;
  while (g_primes[i] <= cap) {
    i++;
  }
  return i;
}

inline void IncrementWheel(int &offset) {
  offset++;
  if (offset >= g_wheelCycleEntries) {
    offset = 0;
  }
}

void MakeWheel(int cap) {
  g_wheelBasePrimes = GetBasePrimes(cap);
  g_wheelCyclePeriod = GetwheelCyclePeriod(cap);
  g_wheelCycleEntries = 0;

  int cur = 0;
  int offset = -1;

  int *wheelBase = 0;

  wheelBase = (int *)malloc(g_wheelCyclePeriod * sizeof(int));
  g_wheelTranslation = (int *)malloc((g_wheelCyclePeriod + 1) * sizeof(int));
  g_wheelOffsets = (int *)malloc(g_wheelCyclePeriod * sizeof(int));
  g_wheelFirstPrime = GetFirstIncludedPrime(cap);

  for (int i = 0; i < g_wheelCyclePeriod; i++) {
    wheelBase[i] = 1;
    for (int j = 2; j <= cap; j++) {
      if (!((i + 1) % j)) {
        wheelBase[i] = 0;
        break;
      }
    }
  }

  while (cur < g_wheelCyclePeriod) {
    if (wheelBase[cur] && cur != 0) {
      g_wheelOffsets[g_wheelCycleEntries] = offset + 1;
      offset = 0;
      g_wheelCycleEntries++;
    } else {
      offset++;
    }
    cur++;
  }

  g_wheelOffsets[g_wheelCycleEntries] = 2;
  g_wheelCycleEntries++;

  int total = 0;

  g_wheelTranslation[0] = 0;
  for (BigInt i = 0; i < g_wheelCyclePeriod; i++) {
    if (i && wheelBase[i - 1]) {
      total++;
    }
    g_wheelTranslation[i + 1] = total;
  }

  free(wheelBase);
}

/* This function calculates how many entries the wheel leaves
in the range from (rangeStart, rangeStop].*/
inline BigInt CountWheelEntries(BigInt rangeStart, BigInt rangeEnd) {
  rangeEnd++;

  int a = rangeStart % g_wheelCyclePeriod;
  int b = rangeEnd % g_wheelCyclePeriod;

  return (rangeEnd - b - rangeStart + a) / g_wheelCyclePeriod *
             g_wheelCycleEntries +
         g_wheelTranslation[b] - g_wheelTranslation[a];
}

void CountHyperbolaLattice_2_Variables(void) {
  BigInt finalBoundary = g_boundary / g_minVarValue;
  BigInt boundaryRoot = (BigInt)(sqrt((double)finalBoundary) + EPSILON);

  /* For if the final two digits happen to be the same. */
  g_latticepoints += g_scale / (g_divisor * (g_divisor + 1));

  /* Leading digit is same, final digit is not. */
  g_latticepoints +=
      (CountWheelEntries(g_minVarValue, finalBoundary / g_minVarValue) - 1) *
      (g_scale / g_divisor);

  /* For if the final two digits happen to be the same,
  but both differ from the previous. */
  g_latticepoints +=
      (CountWheelEntries(g_minVarValue, boundaryRoot) - 1) * (g_scale / 2);

  /* Both digits differ from all other digits - This is the hellish evil
  loop of portentious doom. */

  int curWheelOffset = g_wheelTranslation[g_minVarValue % g_wheelCyclePeriod];

  BigInt curLeadingVar = g_minVarValue + g_wheelOffsets[curWheelOffset];
  BigInt subTotal = 0;

  IncrementWheel(curWheelOffset);

  while (curLeadingVar <= boundaryRoot) {
    subTotal +=
        CountWheelEntries(curLeadingVar, finalBoundary / curLeadingVar) - 1;
    curLeadingVar += g_wheelOffsets[curWheelOffset];
    IncrementWheel(curWheelOffset);
  }

  g_latticepoints += subTotal * g_scale;
}

void CountHyperbolaLattice_3_Variables(BigInt hyperbolaBoundary,
                                       BigInt minVarValue) {
  BigInt maxVarValue = InversePower(hyperbolaBoundary, 3);

  int curWheelOffset = g_wheelTranslation[minVarValue % g_wheelCyclePeriod];

  g_boundary = hyperbolaBoundary;
  g_minVarValue = minVarValue;
  g_scale = g_scaleVal / g_lastScaleDivisor;
  g_divisor = g_lastScaleDivisor + 1;

  CountHyperbolaLattice_2_Variables();

  g_minVarValue += g_wheelOffsets[curWheelOffset];
  IncrementWheel(curWheelOffset);

  g_scale = g_scaleVal;
  g_divisor = 2;
  while (g_minVarValue <= maxVarValue) {
    CountHyperbolaLattice_2_Variables();
    g_minVarValue += g_wheelOffsets[curWheelOffset];
    IncrementWheel(curWheelOffset);
  }
}

void CountHyperbolaLattice_X_Variables(BigInt hyperbolaBoundary,
                                       BigInt minVarValue) {
  BigInt maxVarValue = InversePower(hyperbolaBoundary, g_variablesLeft);

  /* Save global variables that will be restored at end of function */

  BigInt scaleVal = g_scaleVal;
  BigInt lastScaleDivisor = g_lastScaleDivisor;

  int curWheelOffset = g_wheelTranslation[minVarValue % g_wheelCyclePeriod];

  g_variablesLeft--;
  g_lastScaleDivisor = lastScaleDivisor + 1;
  g_scaleVal = scaleVal / lastScaleDivisor;

  if (g_variablesLeft == 3) {
    CountHyperbolaLattice_3_Variables(hyperbolaBoundary / minVarValue,
                                      minVarValue);
  } else {
    CountHyperbolaLattice_X_Variables(hyperbolaBoundary / minVarValue,
                                      minVarValue);
  }

  g_lastScaleDivisor = 2;
  g_scaleVal = scaleVal;
  minVarValue += g_wheelOffsets[curWheelOffset];
  IncrementWheel(curWheelOffset);

  while (minVarValue <= maxVarValue) {
    if (g_variablesLeft == 3) {
      CountHyperbolaLattice_3_Variables(hyperbolaBoundary / minVarValue,
                                        minVarValue);
    } else {
      CountHyperbolaLattice_X_Variables(hyperbolaBoundary / minVarValue,
                                        minVarValue);
    }
    minVarValue += g_wheelOffsets[curWheelOffset];
    IncrementWheel(curWheelOffset);
  }

  /* Restore global variables */

  g_lastScaleDivisor = lastScaleDivisor;
  g_variablesLeft++;
}

BigInt CountHyperbolaLattice(BigInt hyperbolaBoundary, int hyperbolaVariables) {
  g_latticepoints = 0;
  g_variablesLeft = hyperbolaVariables;
  g_lastScaleDivisor = 1;
  g_scaleVal = g_factorial[hyperbolaVariables];

  if (hyperbolaBoundary <
      (BigInt)pow((double)g_wheelFirstPrime, (double)hyperbolaVariables)) {
    return 0;
  }

  switch (hyperbolaVariables) {
  case 1:
    g_latticepoints = CountWheelEntries(g_wheelFirstPrime, hyperbolaBoundary);
    break;
  case 2:
    /* CountHyperbolaLattice_2_Variables expects a number of global variables
    to be initialized when it is called, which generally happens in
    CountHyperbolaLattice_3_Variables.  We have to do it manually here. */
    g_minVarValue = g_wheelFirstPrime;
    g_boundary = g_wheelFirstPrime * hyperbolaBoundary;
    g_scale = 2;
    g_divisor = 1;

    CountHyperbolaLattice_2_Variables();
    break;
  case 3:
    CountHyperbolaLattice_3_Variables(hyperbolaBoundary, g_wheelFirstPrime);
    break;
  default:
    CountHyperbolaLattice_X_Variables(hyperbolaBoundary, g_wheelFirstPrime);
    break;
  }
  return g_latticepoints;
}

BigInt primeCountingFunction(BigInt n) {
  int maxPower = (int)(log((double)n + EPSILON) /
                           log((double)g_wheelFirstPrime + EPSILON) +
                       EPSILON) +
                 1;
  double total = 0.0;

  int oldClock = (int)clock();
  int totalTime = 0;

  for (int curPower = 1; curPower < maxPower; curPower++) {
    if (!g_mu[curPower]) {
      continue;
    }

    BigInt curMax = InversePower(n, curPower);
    double subTotal = 0.0;
    BigInt hyperbolaEntries = 1;
    double sign = 1;

    while (1) {
      double temp = sign / hyperbolaEntries;
      sign *= -1;

      double v2 = (double)CountHyperbolaLattice(curMax, hyperbolaEntries);
      temp *= v2;

      if (temp == 0.0) {
        break;
      }

      subTotal += temp;
      hyperbolaEntries++;

      int newClock = (int)clock();
      totalTime += newClock - oldClock;
      oldClock = newClock;
    }
    subTotal /= curPower * g_mu[curPower];
    total += subTotal;
  }
  total += g_wheelBasePrimes;

  /* the .5 is to prevent truncation errors - but it's clearly sloppy */
  return (BigInt)(total + 0.5);
}

// https://mathoverflow.net/questions/97694/the-relationship-between-the-dirichlet-hyperbola-method-the-prime-counting-func
