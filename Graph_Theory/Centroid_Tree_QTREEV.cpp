#include <bits/stdc++.h>
 
using namespace std;
 
const int maxn=2e5+10;
vector<int> grafo[maxn];//graph
bool rm[maxn];//is removed?
int sub[maxn];//subtree size;

int dfs(int u,int p=-1){//dfs to find sizes;
    sub[u]=1;
    for(auto k:grafo[u]){
        if (k!=p && !rm[k])sub[u]+=dfs(k,u);
    }
    return sub[u];
}

int get_centroid(int u,int size,int p=-1){
    for(auto k:grafo[u]){
        if (k!=p && !rm[k]){
            if (sub[k]*2 > size){//if there is a subtree with size greater than size/2 i can find a centroid there;
                return get_centroid(k,size,u);
            }
        }
    }
    return u;//if there is no such subtree, then i am a centroid;
}

vector<int> centroidtree[maxn];
int pai[maxn];
int root=-1;
int decomposition(int u=0){//doesn't matter which beggining node, calculates the centroid decomposition of the tree;
    int centroid=get_centroid(u,dfs(u));
    if (root==-1)root=centroid;
    rm[centroid]=1;
    for(auto k:grafo[centroid]){
        if (!rm[k]){
            auto l1=decomposition(k);
            pai[l1]=centroid;
            centroidtree[centroid].push_back(l1);
            centroidtree[l1].push_back(centroid);
        }
    }
    return centroid;
}

const int maxlog=17;
int dp[maxn][maxlog];

//quem eh o j ancestral de a em O(logN)
int jmp(int a,int j){
  // dou um bit shift de i, se eu tenho que subir 2^i eu subo 2^i;
  for (int i=0;i<maxlog;i++) if ((j>>i) & 1) {
    a=dp[a][i];
  }
  return a;
}

int cont=0;

//your usual predfs;
int ddepth[maxn];
void predfs(int i=1,int p=-1,int cont=0){
  ddepth[i]=cont;
  for (auto k:grafo[i]) {
    if(k!=p){
      dp[k][0]=i;
      predfs(k,i,cont+1);
    }
  }
}

int lca(int a,int b){

  if (ddepth[a]>ddepth[b]) swap(a,b);
  b=jmp(b,ddepth[b]-ddepth[a]);
  if (a==b) return a;

  for (int i=maxlog-1;i>=0;i--){
    if (dp[a][i]!=dp[b][i]){
      a=dp[a][i];
      b=dp[b][i];
    }
  }
  return dp[a][0];
}
int dist(int a,int b){
  auto lc=lca(a,b);
  return ddepth[a]+ddepth[b]-2*ddepth[lc];
}
int color[maxn];
multiset<int> subarvore[maxn];
int query(int a){
  int i=a;
  int melhor=INT_MAX;
  while(i!=0){
    if(!subarvore[i].empty()){
      auto menor=*subarvore[i].begin();
      melhor=min(melhor,menor+dist(a,i));
    }
    i=pai[i];
  }
  if(melhor==INT_MAX)return -1;
  return melhor;
}
void change(int a){
  if(color[a]==0){
    int i=a;
    while(i!=0){
      //cout<<i<<endl;
      subarvore[i].insert(dist(a,i));
      i=pai[i];
    }
  }
  else{
    int i=a;
    while(i!=0){
      subarvore[i].erase(subarvore[i].find(dist(a,i)));
      i=pai[i];
    } 
  }
  color[a]=1-color[a];
}
signed main(){
  cin.tie(NULL);
  ios_base::sync_with_stdio(false);
    int n;cin>>n;
    for(int i=0;i<n-1;i++){
        int a,b;cin>>a>>b;
        grafo[a].push_back(b);
        grafo[b].push_back(a);
    }
    root=decomposition(1);
    
    //memset(rm,0,sizeof(rm));
    predfs(1);
    //cout<<"a"<<endl;

    for (int j=1;j<maxlog;j++) 
      for (int i=1;i<=n;i++) 
        dp[i][j]=dp[dp[i][j-1]][j-1];
    int q;cin>>q;
    while(q--){
      int a,b;cin>>a>>b;
      //cout<<a<<b<<endl;
      if(a==0){
        change(b);
      }
      else{
        cout<<query(b)<<'\n';
      }
    }
}
