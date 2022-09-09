
struct point {
  flt x, y;
  point(flt x = 0, flt y = 0) : x(x), y(y) {}
  bool operator<(const point &rhs) const {
    return cmp(x, rhs.x) < 0 || (cmp(x, rhs.x) == 0 && cmp(y, rhs.y) < 0);
  }
  bool operator==(const point &rhs) const {
    return cmp(x, rhs.x) == 0 && cmp(y, rhs.y) == 0;
  }
  point operator+(const point &rhs) const {
    return point(x + rhs.x, y + rhs.y);
  }
  point operator-(const point &rhs) const {
    return point(x - rhs.x, y - rhs.y);
  }
  point operator*(const flt k) const { return point(x * k, y * k); }
  point operator/(const flt k) const { return point(x / k, y / k); }
  point operator~() const { // counter clockwise rotate 90 degree
    return point(-y, x);
  }
  flt dot(const point &rhs) const { return x * rhs.x + y * rhs.y; }
  flt det(const point &rhs) const { return x * rhs.y - y * rhs.x; }
  flt norm2() const { return x * x + y * y; }
  flt norm() const { return hypot(x, y); }
  point rot(flt a) const { // counter clockwise rotate A rad
    return point(x * cos(a) - y * sin(a), x * sin(a) + y * cos(a));
  }
  point rot(flt cosa,
            flt sina) const { // counter clockwise rotate using cos/sin
    return point(x * cosa - y * sina, x * sina + y * cosa);
  }
  point trunc(flt a = 1.0) const { return (*this) * (a / this->norm()); }
};
