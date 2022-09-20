BigInt lehmergcd(BigInt p, BigInt q) {
  int cmpres = cmp(p, q);
  if (cmpres == 0)
    return p;
  if (cmpres < 0)
    swap(p, q);
  unsigned long long x, y, z, num1, den1, w1, num2, den2, w2, e, f, xn, yn, t;
  unsigned int a, b, c, d, w;
  bool needlongdiv;
  int parity;
  int nlong = 0, nlehmer = 0, clehmer, nit = 0;
  while (true) {
    if (SZ(q.d) == 0)
      return p;
    else if (SZ(p.d) <= 2)
      break;
    else
      needlongdiv = false;
    if (SZ(p.d) - SZ(q.d) >= 2)
      needlongdiv = true;
    if (!needlongdiv) {
      extractleadingbits(p, q, x, y);
      if (y <= BIGINTMASK || x == y)
        needlongdiv = true;
      if (x == ((((unsigned long long)BIGINTMASK) << BIGINTBITS) | BIGINTMASK))
        x >>= 1, y >>= 1;
    }
    if (!needlongdiv) {
      num1 = x, den1 = y + 1, num2 = x + 1, den2 = y, w1 = num1 / den1,
      w2 = num2 / den2;
      if (w1 != w2 || w1 > BIGINTMASK)
        needlongdiv = true;
      else
        w = w1;
    }
    if (!needlongdiv) {
      a = 0, b = 1, c = 1, d = w, z = x - w * y, x = y, y = z, parity = 0,
      clehmer = 1;
      while (true) {
        if (parity == 0) {
          if (y == d)
            break;
          num1 = x - a, den1 = y + c, num2 = x + b, den2 = y - d;
        }
        if (parity == 1) {
          if (y == c)
            break;
          num1 = x - b, den1 = y + d, num2 = x + a, den2 = y - c;
        }
        w1 = num1 / den1, w2 = num2 / den2;
        if (w1 != w2 || w1 > BIGINTMASK)
          break;
        else
          w = w1;
        e = a + w * c, f = b + w * d, z = x - w * y;
        if (e > BIGINTMASK || f > BIGINTMASK)
          break;
        else
          a = c, c = e, b = d, d = f, x = y, y = z, parity = 1 - parity,
          ++clehmer;
      }
    }
    if (!needlongdiv && b != 0) {
      // printf("lehmer step (%s,%s) -> %u %u %u %u
      // (parity=%d)\n",format(p).c_str(),format(q).c_str(),a,b,c,d,parity);
      x = 0, y = 0, xn = 0, yn = 0, nlehmer += clehmer, ++nit;
      while (SZ(q.d) < SZ(p.d))
        q.d.PB(0);
      for (int i = 0; i < SZ(p.d); ++i) {
        unsigned long long cp = p.d[i], cq = q.d[i];
        if (parity == 0)
          x += cq * b, xn += cp * a, y += cp * c, yn += cq * d;
        else
          x += cp * a, xn += cq * b, y += cq * d, yn += cp * c;
        t = min(x, xn), x -= t, xn -= t, t = min(y, yn), y -= t, yn -= t;
        if (xn == 0)
          p.d[i] = x & BIGINTMASK, x >>= BIGINTBITS;
        else if ((xn & BIGINTMASK) == 0)
          p.d[i] = 0, xn >>= BIGINTBITS;
        else
          p.d[i] = BIGINTMASK - (xn & BIGINTMASK) + 1, xn >>= BIGINTBITS, ++xn;
        if (yn == 0)
          q.d[i] = y & BIGINTMASK, y >>= BIGINTBITS;
        else if ((yn & BIGINTMASK) == 0)
          q.d[i] = 0, yn >>= BIGINTBITS;
        else
          q.d[i] = BIGINTMASK - (yn & BIGINTMASK) + 1, yn >>= BIGINTBITS, ++yn;
        // if(parity==0) x+=cq*b-cp*a,y+=cp*c-cq*d; else
        // x+=cp*a-cq*b,y+=cq*d-cp*c;
        // p.d[i]=x&BIGINTMASK,x>>=BIGINTBITS,q.d[i]=y&BIGINTMASK,y>>=BIGINTBITS;
      }
      // printf("=> (%s,%s)
      // [%llx,%llx,%llx,%llx]\n",format(p).c_str(),format(q).c_str(),x,y,xn,yn);
      assert(x == 0 && y == 0 && xn == 0 && yn == 0);
      normalize(p);
      normalize(q);
    } else {
      BigInt r = p % q;
      p = q, q = r;
      ++nlong, ++nit;
    }
  }
  x = (((unsigned long long)(SZ(p.d) == 2 ? p.d[1] : 0)) << BIGINTBITS) |
      p.d[0];
  y = (((unsigned long long)(SZ(q.d) == 2 ? q.d[1] : 0)) << BIGINTBITS) |
      q.d[0];
  while (y != 0) {
    z = x % y, x = y, y = z;
  }
  // printf("nlong=%d nlehmer=%d nit=%d =>
  // %.2lf\n",nlong,nlehmer,nit,1.0*nlehmer/(nit-nlong));
  return BigInt(x);
}
