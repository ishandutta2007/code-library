#include <bits/stdc++.h>
using namespace std;
#define LL long long
#define EPS 0.00000001
#define PI 2*acos(0.0)
#define Fout freopen("output.txt","w",stdout)
#define Fin freopen("input.txt","r",stdin)
#define MS(ARRAY,VALUE) memset(ARRAY,VALUE,sizeof(ARRAY))
#define RT printf("Run Time : %0.3lf seconds\n", clock()/(CLOCKS_PER_SEC*1.0))
LL prime[100010];
bool mark[1000030];
LL nprime=0;

template <class X> void sieve(X n)
{
    X i,j,limit = sqrt(n)+2;

    mark[1] = 1;

    for(i=4;i<=n;i+=2) mark[i]=1;
    prime[++nprime]=2;
    for(i=3;i<=n;i++)
    {
        if(!mark[i])
        {
            prime[++nprime]=i;

            if(i<=limit)
            {
                for(j=i*i;j<=n;j+=i*2)
                {
                    mark[j] = 1;
                }
            }
        }
    }
}
template <class X> X NOD(X n)
{
    if(n==1) return 1;
    X s = sqrt(n)+2,res=1,cnt=0;
    X i = 1;
    while(prime[i]<=s)
    {
        if(n%prime[i]==0)
        {
            n=n/prime[i];
            cnt++;
            if(n==1) {res=res*(cnt+1);break;}
        }
        else
        {
            i++;
            res=res*(cnt+1);
            cnt=0;
        }
    }
    if(n != 1) res=res*2;
    return res;
}

int main()
{
    //Fout;
    sieve(1000020);
    //int n;
    //for(int i=2;i<=10000;i++) cout <<i<<" -> "<< NOD(i) <<endl;
    cout <<NOD(1000000000000LL)<<endl;
    //RT;
    return 0;
}
