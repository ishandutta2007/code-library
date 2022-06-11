#include <bits/stdc++.h>
using namespace std;
#define LL long long
#define INF INT_MAX
#define EPS 0.00000001
#define PI 2*acos(0.0)
#define MOD 1000000007
#define ck(XX) cout<<XX<<endl
#define Set(XX,POS) XX|(1<<POS)
#define reset(XX,POS) XX&(~(1<<POS))
#define check(XX,POS) (bool)(XX&(1<<POS))
#define toggle(XX,POS) (XX^(1<<POS))
#define Fin freopen("input.txt","r",stdin)
#define Fout freopen("output.txt","w",stdout)
#define valid(X,Y,R,C) X>=0 && X<R && Y>=0 && Y<C
#define MS(ARRAY,VALUE) memset(ARRAY,VALUE,sizeof(ARRAY))
#define RT printf("Run Time : %0.3lf seconds\n", clock()/(CLOCKS_PER_SEC*1.0))
int prime[10000000];
int mark[100000100/32];
int nprime=0;
void Bsieve(int n)
{
	int i,j,limit = sqrt(n)+2;

	mark[1/32]=Set(mark[1/32], 1%32);

	for(i=4; i<=n; i+=2) mark[i/32]=Set(mark[i/32], i%32);
	prime[++nprime] = 2;

	for(i=3; i<=n; i+=2)
	{
		if(check(mark[i/32],i%32)==0)
		{
			prime[++nprime] = i;
			if(i <= limit)
			{
				for(j=i*i; j<=n; j+=2*i)
				{
					mark[j/32]=Set(mark[j/32], j%32);
				}
			}
		}
	}
	return;
}




int main()
{
	Bsieve(100000100);
	ck(nprime);
	return 0;
}
