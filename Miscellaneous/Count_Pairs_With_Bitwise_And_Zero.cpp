// https://codeforces.com/blog/entry/20174
// O(A_i_max log A_i_max)

#include <cstdio>
    #include <vector>
    using namespace std;
     
    const int MAXN = 1048576;
    long long solve(vector<int> &contain)
    {
    	long long ans = 0;
    	for (int j = 0; (1<<j) < MAXN; j++) {
    		for (int i = 0; i < MAXN; i++) if (!(i & (1<<j))) {
    			contain[i] += contain[i^(1<<j)];
    		}
    	}
    	for (int i = 0; i < MAXN; i++) {
    		if (__builtin_popcount(i) & 1) ans -= (contain[i]*1ll*(contain[i]-1));
    		else ans += (contain[i]*1ll*(contain[i]-1));
    	}
    	return ans;
    }
    void tmain()
    {
    	int N;
    	vector<int> contain(MAXN);
    	scanf("%d",&N);
    	for(int i=0;i<N;i++)
    	{
    		int t;
    		scanf("%d",&t);
    		contain[t]++;
    	}

    	printf("%lld\n",contain[0]+solve(contain));
    }
    int main()
    {
    	int T;
    	scanf("%d",&T);
    	while(T--)
    		tmain();
    	return 0;
    }

