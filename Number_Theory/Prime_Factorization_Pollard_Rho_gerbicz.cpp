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
  long long int a, n, G, szorzat, tomb[293]; // nagy pr�mvari�nshoz kell a t�mb

  // 2^16-ig a pr�mek gener�l�sa Erathoszteneszi szit�val
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

  printf("Pollard p-1 faktorizacios modszere legfeljebb 63 bites szamokra!\n");
  printf("Nagy primvarianssal.\n");
  printf("Lnkot csak az elso/masodik lepes utan szamol a program!\n");
  printf("p (prim)osztot talal a Pollard, ha p-1=Q*R, ahol Q minden\n");
  printf("primhatvanyosztoja <=B1, tovabba R=1 vagy R prim es (B1<) R<=B2\n");
  printf("Ha lnko=n lenne, akkor probalj kisebb B1/B2 ertekkel!\n");
  printf("B2<=B1 eseten masodik lepest mar nem szamol a program.\n");
  printf("Korlatok: 2<=B1,B2<2^31\n\n");
  printf("Kerem a szamot: (0-ra kilep a program)\n");

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

    a = 2; // az alap a Pollard m�dszerben, nem l�nyeges
    for (I = 0; I <= B1; I += 65536) {
      felso = I + 65536;
      if (felso > B1)
        felso = B1 + 1;
      L = felso - I;
      // szita a k�vetkez� pr�mek megtal�l�s�ra
      for (i = 0; i < L; i++)
        isprime[i] = 1;
      if (I == 0)
        isprime[0] = 0, isprime[1] = 0;
      for (i = 0; (i < 6542) && (prime[i] * prime[i] <= felso); i++) {
        p = prime[i];
        elso = ((I + p - 1) / p) * p;
        if (elso <= p)
          elso = 2 * p; // kis primeket ne huzzuk ki
        for (j = elso - I; j < L; j += p)
          isprime[j] = 0;
      }
      db = 0;
      for (i = 0; i < L; i++)
        if (isprime[i]) {
          p = I + i;
          q = p;
          while (q <= B1 / p)
            q *= p; // B1-ig hatvanyozzuk p-t
          pr[db] = q;
          db++;
        }
      // a Pollard m�dszer lelke
      for (i = 0; i < db; i++)
        a = powmod(a, pr[i], n);
    }

    G = gcd(a - 1, n);
    if ((G != 1) && (G != n)) {
      printf("Nemtrivialis osztot talalt a p-1 modszer az elso lepesben:\n");
      printf("n osztoja: %lld\n", G);
      tovabb = 0;
    } else if (G == n) {
      printf("lnko=n, probalj kisebb B1-el!\n");
      tovabb = 0;
    } else {
      printf("lnko=1 az elso lepesben.\n");
      tovabb = 1;
    }
    if (tovabb && (B2 > B1)) {

      printf("Nagy prim varians inditasa.\n");

      // tomb[i]=a^i mod n
      tomb[0] = 1;
      for (i = 1; i <= 292;
           i++) // maxim�lis gap=292 a legfeljebb 31 bites pr�mekre
        tomb[i] = mulmod(a, tomb[i - 1], n);
      a = powmod(a, B1, n);
      elozo = B1;
      szorzat = 1;

      for (I = B1 + 1; I <= B2; I += 65536) {
        felso = I + 65536;
        if (felso > B2)
          felso = B2 + 1;
        L = felso - I;
        // szita a k�vetkez� pr�mek megtal�l�s�ra
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

        // Pollard m�dszer m�sodik l�p�se
        for (i = 0; i < db; i++) {
          a = mulmod(a, tomb[pr[i] - elozo], n);
          szorzat = mulmod(szorzat, a - 1, n);
          elozo = pr[i];
        }
      }
      G = gcd(szorzat, n);
      if ((G != 1) && (G != n)) {
        printf(
            "Nemtrivialis osztot talalt a p-1 modszer a masodik lepesben:\n");
        printf("n osztoja: %lld\n", G);
      } else if (G == n) {
        printf("lnko=n, probalj kisebb B2-el!\n");
      } else {
        printf("lnko=1 a masodik lepesben is.\n");
      }
    }
    printf("\n");
    printf("Kerem a szamot: (0-ra kilep a program)\n");
  }
  return 0;
}
