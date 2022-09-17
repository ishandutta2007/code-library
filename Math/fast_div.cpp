#define DIVIDE(n, d) ((FastInt)((double)(n)/d))
// or
auto fast_div=[](const ll& a,const int& b)->ll {return double(a)/b+1e-9;};
