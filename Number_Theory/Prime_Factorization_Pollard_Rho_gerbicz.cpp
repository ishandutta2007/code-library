// credits: Robert Gerbicz
#include <assert.h>
#include <stdio.h>

long long int mulmod(long long int a, long long int b, long long int p)
// (a*b)%p kisz�m�t�sa, klasszikus
{
  long long int y = (long long int)((double)a * (double)b / p + 0.5);
  long long int r = a * b - y * p;
  r %= p;
  if (r < 0)
    r += p;

  return r;
}

long long int powmod(long long int alap, long long int kitevo,
                     long long int modulus) {
  // alap^kitevo mod modulus-sal t�r vissza, itt kitevo>=0
  long long int eredmeny = 1, hatvany = alap;

  while (kitevo) {
    if (kitevo & 1)
      eredmeny = mulmod(eredmeny, hatvany, modulus);
    hatvany = mulmod(hatvany, hatvany, modulus);
    kitevo >>= 1;
  }
  return eredmeny;
}

long long int gcd(long long int a,
                  long long int b) { // lnko(a,b)-vel t�r vissza
  // az oszt�sok nagy r�sz�t elker�l� gyors algoritmus, b�r ezt
  // most legfeljebb 2-szer sz�molunk �gy a sebess�ge nem l�nyeges
  if (a < 0)
    a = -a;
  if (b < 0)
    b = -b;
  if (a == 0)
    return b;
  if (b == 0)
    return a;

  long long int c;

  while (b > 0) {
    if (a >= b) {
      a -= b;
      if (a >= b) {
        a -= b;
        if (a >= b) {
          a -= b;
          if (a >= b) {
            a -= b;
            if (a >= b) {
              a -= b;
              if (a >= b) {
                a -= b;
                if (a >= b) {
                  a -= b;
                  if (a >= b) {
                    a -= b;
                    if (a >= b)
                      a %= b;
                  }
                }
              }
            }
          }
        }
      }
    }
    c = a, a = b, b = c;
  }
  return a;
}

int main() {

  unsigned int db, tovabb, B1, B2, I, i, j, L, p, q, elso, felso, elozo,
      isprime[65536], prime[6542], pr[6542]; // pi(2^16)=6542
  long long int a, n, G, szorzat, tomb[293]; // for large prime varins, you need the array
  // Generation of prims up to 2^16 with the Sieve of Eratosthenes
  for (i = 0; i < 65536; i++)
    isprime[i] = (i > 1);
  for (i = 0; i < 256; i++) {
    if (isprime[i]) {
      for (j = i * i; j < 65536; j += i)
        isprime[j] = 0;
    }
  }
  db = 0;
  for (i = 0; i < 65536; i++)
    if (isprime[i])
      prime[db] = i, db++;

  printf("Pollard's p-1 factorization mod for up to 63-bit numbers!\n");
  printf("With large prime variance.\n");
  printf("The program only collects a link after the first/second page!\n");
  printf("Pollard finds p (prim) division if p-1=Q*R, where Q is all\n");
  printf("prime sixty divisor <=B1, then R=1 or R prime and (B1<) R<=B2\n");
  printf("If lnko=n, then try smaller B1/B2 values!\n");
  printf("In the case of B2<=B1, the program does not generate a second error.\n");
  printf("Constraints: 2<=B1,B2<2^31\n\n");
  printf("Please enter the number: (program exits to 0)\n");

  while (scanf("%lld", &n) != EOF) {
    if (n == 0)
      return 0;

    if (n < 0)
      n = -n;

    printf("B1 korlat: ");
    scanf("%u", &B1);
    printf("B2 korlat: ");
    scanf("%u", &B2);
    assert(B1 <= 2147483647 && B2 <= 2147483647);

    a = 2; // the base in the Pollard method is not essential
    for (I = 0; I <= B1; I += 65536) {
      felso = I + 65536;
      if (felso > B1)
        felso = B1 + 1;
      L = felso - I;
      // sieve to find the next items
      for (i = 0; i < L; i++)
        isprime[i] = 1;
      if (I == 0)
        isprime[0] = 0, isprime[1] = 0;
      for (i = 0; (i < 6542) && (prime[i] * prime[i] <= felso); i++) {
        p = prime[i];
        elso = ((I + p - 1) / p) * p;
        if (elso <= p)
          elso = 2 * p; // do not draw small primes
        for (j = elso - I; j < L; j += p)
          isprime[j] = 0;
      }
      db = 0;
      for (i = 0; i < L; i++)
        if (isprime[i]) {
          p = I + i;
          q = p;
          while (q <= B1 / p)
            q *= p; // multiply p to B1
          pr[db] = q;
          db++;
        }
      // the soul of Pollard medicine
      for (i = 0; i < db; i++)
        a = powmod(a, pr[i], n);
    }

    G = gcd(a - 1, n);
    if ((G != 1) && (G != n)) {
      printf("The p-1 method found a non-trivial division in the first sheet:\n");
      printf("divisor of n: %lld\n", G);
      tovabb = 0;
    } else if (G == n) {
      printf("lnko=n, try smaller B1!\n");
      tovabb = 0;
    } else {
      printf("lnko=1 in the first sheet.\n");
      tovabb = 1;
    }
    if (tovabb && (B2 > B1)) {

      printf("High prime variance index.\n");

      // tomb[i]=a^i mod n
      tomb[0] = 1;
      for (i = 1; i <= 292; i++) // maximum gap=292 for up to 31-bit pr�ms
        tomb[i] = mulmod(a, tomb[i - 1], n);
      a = powmod(a, B1, n);
      elozo = B1;
      szorzat = 1;

      for (I = B1 + 1; I <= B2; I += 65536) {
        felso = I + 65536;
        if (felso > B2)
          felso = B2 + 1;
        L = felso - I;
        // sieve to find the next items
        for (i = 0; i < L; i++)
          isprime[i] = 1;
        for (i = 0; (i < 6542) && (prime[i] * prime[i] <= felso); i++) {
          p = prime[i];
          elso = ((I + p - 1) / p) * p;
          if (elso <= p)
            elso = 2 * p;
          for (j = elso - I; j < L; j += p)
            isprime[j] = 0;
        }
        db = 0;
        for (i = 0; i < L; i++)
          if (isprime[i])
            pr[db] = I + i, db++;

        // Pollard's second step
        for (i = 0; i < db; i++) {
          a = mulmod(a, tomb[pr[i] - elozo], n);
          szorzat = mulmod(szorzat, a - 1, n);
          elozo = pr[i];
        }
      }
      G = gcd(szorzat, n);
      if ((G != 1) && (G != n)) {
        printf("The p-1 mod scored a non-trivial split in the second round:\n");
        printf("n refuge: %lld\n", G);
      } else if (G == n) {
        printf("lnko=n, try a smaller B2!\n");
      } else {
        printf("lnko=1 also in the second lepes.\n");
      }
    }
    printf("\n");
    printf("Kerem and chamotte: (0-ra kilep and program)\n");
  }
  return 0;
}
