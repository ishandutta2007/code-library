
pair <int, int> reduce_fraction(int num, int den){
  int g = gcd(num, den);
  return make_pair(num/g,den/g);
}