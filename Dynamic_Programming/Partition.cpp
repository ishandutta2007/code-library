#include <iostream>
const int MOD = 1000000;
int* p = new int [100000+7];

int partition(int n)
{
	int sum = 0;
	int k = 4;
	int a = 2;
	int b = 1;
	int sgn = 1;
	
	while (n >= a)
	{
		sum += sgn*(p[n-a] + p[n-b]);
		sum = sum;
		a += k + 1;
		b += k;
		sgn *= -1; 
		k+=3;
	}
	
	if (n >= b)
	{
		sum += sgn*p[n-b];
	}
	return sum % MOD;
}

int main()
{
	p[0] = 1;
	p[1] = 1;
	int n = 1;
	for (int n = 2; n <= 10000; i++){
		n++;
		p[n] = partition(n);
	}
	std::cout << p[10000] << "\n";
	return 1;
}
