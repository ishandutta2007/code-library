
long double non_modular_power(double a, int b) {
  if (b == 0)
    return 1.;
  if (b & 1)
    return non_modular_power(a, b - 1) * a;
  long double half = non_modular_power(a, b >> 1);
  return half * half;
}
