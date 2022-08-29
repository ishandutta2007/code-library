// See codechef FLIPCOIN for an example
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>

/*
Output of segment tree range queries are modelled as elements of a monoid.
Updates on the segment tree are functions on the monoid.

A Monoid is a class implementing the following methods:
* Monoid();  // constructor: identity element of monoid
* Monoid(const Monoid&, const Monoid&);  // constructor: element created by
combining 2 elements
* void print(FILE* fp);

A Function is a class implementing the following methods:
* Function();  // constructor: identity function
* Function(const Function& f1, const Function& f2);  // constructor: function
composition (f1.f2)
* bool is_identity() const;  // whether this function is the identity function
* Monoid operator()(const Monoid&);  // function definition
* void print(FILE* fp);
*/

int segtree_size(int n) {
  // 2^(ceil(lg(n)) + 1) - 1
  n -= 1;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return 2 * n + 1;
}

template <class M, class F> class SegTree {
  // M is a monoid
  // F is a function
public:
  typedef M value_type;
  typedef F func_type;

  static void identity_check() {
    if (!F().is_identity()) {
      fprintf(stderr,
              "SegTree: function's default constructor is not identity\n");
      std::exit(1);
    }
  }

private:
  int n;
  std::vector<M> values;
  std::vector<F> pends;

public:
  explicit SegTree(int _n)
      : n(_n), values(segtree_size(n)), pends(values.size()) {
    identity_check();
  }
  explicit SegTree(int _n, const M &x)
      : n(_n), values(segtree_size(n)), pends(values.size()) {
    identity_check();
    assign(x);
  }
  explicit SegTree(int _n, const M *a)
      : n(_n), values(segtree_size(n)), pends(values.size()) {
    identity_check();
    assign(a);
  }
  explicit SegTree(const std::vector<M> &v)
      : n(v.size()), values(segtree_size(v.size())), pends(values.size()) {
    identity_check();
    assign(v);
  }

private:
  void assign_values(int root, int first, int last, const M *a) {
    // root has the node number, first and last have the array indices.
    if (first == last) {
      values[root] = a[first];
    } else {
      int left = 2 * root + 1;
      int mid = (first + last) / 2;
      assign_values(left, first, mid, a);
      assign_values(left + 1, mid + 1, last, a);
      values[root] = M(values[left], values[left + 1]);
    }
  }

  void assign_values(int root, int first, int last, const M &x) {
    // root has the node number, first and last have the array indices.
    if (first == last) {
      values[root] = x;
    } else {
      int left = 2 * root + 1;
      int mid = (first + last) / 2;
      assign_values(left, first, mid, x);
      assign_values(left + 1, mid + 1, last, x);
      values[root] = M(values[left], values[left + 1]);
    }
  }

public:
  void assign(const M *a) {
    pends.assign(values.size(), F());
    assign_values(0, 0, n - 1, a);
  }

  void assign(const M &x) {
    pends.assign(values.size(), F());
    assign_values(0, 0, n - 1, x);
  }

  void assign(const std::vector<M> &v) {
    if (v.size() < n) {
      fprintf(stderr, "SegTree: vector input to assign is too short\n");
      std::exit(3);
    }
    assign(v.data);
  }

private:
  void print(int root, int first, int last, FILE *fp, int level) const {
    for (int i = 0; i < level; ++i) {
      fprintf(fp, "  ");
    }
    fprintf(fp, "%d[%d-%d]: ", root, first, last);
    values[root].print(fp);
    if (!pends[root].is_identity()) {
      fprintf(fp, ": ");
      pends[root].print(fp);
    }
    fprintf(fp, "\n");
    if (first != last) {
      int left = 2 * root + 1;
      int mid = (first + last) / 2;
      print(left, first, mid, fp, level + 1);
      print(left + 1, mid + 1, last, fp, level + 1);
    }
  }

public:
  void print(FILE *fp) const { print(0, 0, n - 1, fp, 0); }

private:
  void propagate(int root, int first, int last) {
    if (!pends[root].is_identity()) {
      values[root] = pends[root](values[root]);
      if (first != last) {
        int left = 2 * root + 1;
        pends[left] = F(pends[root], pends[left]);
        pends[left + 1] = F(pends[root], pends[left + 1]);
      }
      pends[root] = F();
    }
  }

  M query(int root, int first, int last, int i, int j) {
    if (i > last || j < first) {
      return M();
    }
    propagate(root, first, last);
    if (i <= first && last <= j) {
      return values[root];
    } else {
      int left = 2 * root + 1;
      int mid = (first + last) / 2;
      return M(query(left, first, mid, i, j),
               query(left + 1, mid + 1, last, i, j));
    }
  }

public:
  M query(int i, int j) { return query(0, 0, n - 1, i, j); }

private:
  void update(int root, int first, int last, int i, int j, const F &f) {
    if (i > last || j < first) {
      propagate(root, first, last);
    } else if (i <= first && last <= j) {
      pends[root] = F(f, pends[root]);
      propagate(root, first, last);
    } else {
      propagate(root, first, last);
      int left = 2 * root + 1;
      int mid = (first + last) / 2;
      update(left, first, mid, i, j, f);
      update(left + 1, mid + 1, last, i, j, f);
      values[root] = M(values[left], values[left + 1]);
    }
  }

public:
  void update(int i, int j, const F &f) { update(0, 0, n - 1, i, j, f); }
};
