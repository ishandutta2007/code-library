#include<bits/stdc++.h>
using namespace std;

const int N = 3e5 + 9;

int a[N];
struct ST {
  int t[4 * N];
  static const int inf = 1e9;
  ST() {
    memset(t, 0, sizeof t);
  }
  void build(int n, int b, int e) {
    if (b == e) {
      t[n] = a[b];
      return;
    }
    int mid = (b + e) >> 1, l = n << 1, r = l | 1;
    build(l, b, mid);
    build(r, mid + 1, e);
    t[n] = max(t[l], t[r]);
  }
  
  //####################################  NEW ####################################
  void buildNew(int n,int b,int e){
    if(b==e){
      t[b]=a[e]; return;
    }
    int cl=2*n,cr=2*n+1,m=(b+e)/2;
    build(cl,b,m);
    build(cr,m+1,e);
    t[n]=max(t[cl],t[cr]);
  }
  //####################################  NEW ####################################
  
  
  
  
  

  
  void update(int n, int b, int e, int i, int x) {// update at index i with x
    if (b > i || e < i) return;
    if (b == e && b == i) {
      t[n] = x;
      return;
    }
    int mid = (b + e) >> 1, l = n << 1, r = l | 1;
    upd(l, b, mid, i, x);
    upd(r, mid + 1, e, i, x);
    t[n] = max(t[l], t[r]);
  }
  //####################################  NEW ####################################
  void updNew(int n,int b,int e,int i,int x){
    if(b==e){
      t[n]=x; return;
    }
    int cl=2*n,cr=2*n+1,m=(b+e)/2;
    if(i<=m) upd(cl,b,m,i,x);
    else upd(cr,m+1,e,i,x);
    t[n]=max(t[cl],t[cr]);
  }
  //####################################  NEW ####################################
  
  
  
  
  
  
  
  int query(int n, int b, int e, int i, int j) {    // [i,j] range in which we want to qry
    if (j<b || e < i) return -inf;
    if (i<=b && e <= j) return t[n];
    int mid = (b + e) >> 1, l = n << 1, r = l | 1;
    int L = query(l, b, mid, i, j);
    int R = query(r, mid + 1, e, i, j);
    return max(L, R);
  }
  //####################################  NEW ####################################
  int qryNew(int n,int b,int e,int i,int j){
    if(b==i&&e==j) return t[n];
    int cl=2*n,cr=2*n+1,m=(b+e)/2;
    if(j<=m) return qry(cl,b,m,i,j);
    else if(i>m)return qry(cr,m+1,e,i,j);
    else{
        int L=qry(cl,b,m,i,m),R=qry(cr,m+1,e,m+1,j);
        return max(L,R);
    }
  }
  //####################################  NEW ####################################
}t;






int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  // qry/upd..(1,1,n,i,j)
  return 0;
}
