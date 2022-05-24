
static int jacobi_iu(IV in, UV m) {
  int j = 1;
  UV n = (in < 0) ? -in : in;

  if (m <= 0 || (m % 2) == 0)
    return 0;
  if (in < 0 && (m % 4) == 3)
    j = -j;
  while (n != 0) {
    while ((n % 2) == 0) {
      n >>= 1;
      if ((m % 8) == 3 || (m % 8) == 5)
        j = -j;
    }
    {
      UV t = n;
      n = m;
      m = t;
    }
    if ((n % 4) == 3 && (m % 4) == 3)
      j = -j;
    n = n % m;
  }
  return (m == 1) ? j : 0;
}