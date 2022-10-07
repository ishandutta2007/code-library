//package ecr012;https://codeforces.com/problemset/problem/665/F

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.InputMismatchException;

public class F3 {
	InputStream is;
	PrintWriter out;
	String INPUT = "";
	
	public static int[] sieveEratosthenes(int n) {
		if (n <= 32) {
			int[] primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
			for (int i = 0; i < primes.length; i++) {
				if (n < primes[i]) {
					return Arrays.copyOf(primes, i);
				}
			}
			return primes;
		}

		int u = n + 32;
		double lu = Math.log(u);
		int[] ret = new int[(int) (u / lu + u / lu / lu * 1.5)];
		ret[0] = 2;
		int pos = 1;

		int[] isnp = new int[(n + 1) / 32 / 2 + 1];
		int sup = (n + 1) / 32 / 2 + 1;

		int[] tprimes = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
		for (int tp : tprimes) {
			ret[pos++] = tp;
			int[] ptn = new int[tp];
			for (int i = (tp - 3) / 2; i < tp << 5; i += tp)
				ptn[i >> 5] |= 1 << (i & 31);
			for (int j = 0; j < sup; j += tp) {
				for (int i = 0; i < tp && i + j < sup; i++) {
					isnp[j + i] |= ptn[i];
				}
			}
		}

		// 3,5,7
		// 2x+3=n
		int[] magic = { 0, 1, 23, 2, 29, 24, 19, 3, 30, 27, 25, 11, 20, 8, 4, 13, 31, 22, 28, 18, 26, 10, 7, 12, 21, 17,
				9, 6, 16, 5, 15, 14 };
		int h = n / 2;
		for (int i = 0; i < sup; i++) {
			for (int j = ~isnp[i]; j != 0; j &= j - 1) {
				int pp = i << 5 | magic[(j & -j) * 0x076be629 >>> 27];
				int p = 2 * pp + 3;
				if (p > n)
					break;
				ret[pos++] = p;
				if ((long) p * p > n)
					continue;
				for (int q = (p * p - 3) / 2; q <= h; q += p)
					isnp[q >> 5] |= 1 << q;
			}
		}

		return Arrays.copyOf(ret, pos);
	}

	int[] primes = sieveEratosthenes(400000);
	
	void solve()
	{
		long n = nl();
		long ans = 0;
		if(n <= 10) {
			if(n >= 6)ans++;
			if(n >= 8)ans++;
			if(n >= 10)ans++;
			out.println(ans);
			return;
		}
		
		MeisselLehmer2 ml = new MeisselLehmer2(n/2);
		for(int i = 0;i < primes.length;i++) {
			long plus = ml.pi(n/primes[i]) - (i+1);
			if(plus <= 0)break;
			ans += plus;
		}
		
		for(int p : primes) {
			if((long)p*p*p <= n) {
				ans++;
			}else {
				break;
			}
		}
		out.println(ans);
	}
	
	public static class MeisselLehmer2 {
		int[] primes;
//		int[] cump;
		int[] xcums;
		long[] isp;
		
		int cachelimit;
//		Map<Long, Long> picache;
//		Map<Long, Long> phicache;
		LongHashCounterL picache;
		LongHashCounterL phicache;
		boolean cachePhi = false;
		
		int[][] cumps;
//		int[] FP = {2,3,5,7,11,13,17};
		int[] FP = {2,3,5,7,11,13,17,19};
		int M = FP.length;
		
		int B = 100;
		// 10^14 -> 200
		// 10^14/8 -> 100~200

		public MeisselLehmer2(long n)
		{
			assert n >= 4;
//			picache = new HashMap<>();
//			phicache = new HashMap<>();
			picache = new LongHashCounterL();
			phicache = new LongHashCounterL();
			
			int s = (int)Math.sqrt(n)*B; // not overflow!
			cachelimit = s;
			this.primes = sieveEratosthenes(s);
			xcums = new int[(s>>>6)+1];
			isp = new long[(s>>>6)+1];
			for(int p : primes)isp[p>>>6] |= 1L<<p;
			for(int i = 0;i < xcums.length-1;i++)xcums[i+1] = xcums[i] + Long.bitCount(isp[i]);
			
			cumps = new int[M][];
			boolean[] b = {true};
			int len = 1;
			for(int i = 0;i < M;i++){
				len *= FP[i];
				boolean[] c = new boolean[len];
				for(int j = 0;j < FP[i];j++){
					System.arraycopy(b, 0, c, j*b.length, b.length);
				}
				for(int k = 0;k < len;k+=FP[i])c[k] = false;
				cumps[i] = new int[len+1];
				for(int j = 0;j < len;j++)cumps[i][j+1] = cumps[i][j] + (c[j] ? 1 : 0);
				FP[i] = len;
				b = c;
			}
		}
		
		private long sqrt(long n)
		{
			long i = Math.max(0, (long)Math.sqrt(n)-2);
			while(i*i <= n)i++;
			return i-1;
		}
		
		private long cbrt(long n)
		{
			long i = Math.max(0, (long)Math.cbrt(n)-2);
			while(i*i*i <= n)i++;
			return i-1;
		}
		
		public long pi(long x)
		{
			if(x <= cachelimit){
				int ix = (int)x;
				return xcums[ix>>>6] + Long.bitCount(isp[ix>>>6]<<~ix);
//				return cump[(int)x];
			}
			if(picache.containsKey(x))return picache.get(x);
			
			int A = (int)pi(sqrt(sqrt(x)));
			long ret = A + phi(x, A) - P2(x, A) - P3(x, A) - 1;
			picache.put(x, ret);
			return ret;
		}
		
		private long phi(long x, int A)
		{
			if(A > 0 && A-1 < M){
				return cumps[A-1][FP[A-1]] * (x/FP[A-1]) + cumps[A-1][(int)(x%FP[A-1]+1)];
			}
			if(A > 0 && x <= (long)primes[A-1]*primes[A-1]){
				return pi(x) - A + 1;
			}
//			long code = x<<13|A;
//			if(phicache.containsKey(code))return phicache.get(code);
			
			long ret = x;
			for(int i = A-1;i >= 0;i--)ret -= phi(x/primes[i], i);
			
//			if(cachePhi)phicache.put(code, ret);
			return ret;
		}
		
		private long P2(long x, int A)
		{
			int B = (int)pi(sqrt(x));
			long ret = 0;
			for(int i = A;i < B;i++){
				ret += pi(x/primes[i]);
			}
			ret -= (long)(B-A)*(B+A-1)/2;
			return ret;
		}
		
		private long P3(long x, int A)
		{
			int C = (int)pi(cbrt(x));
			long ret = 0;
			for(int i = A;i < C;i++){
				long xi = x/primes[i];
				int B = (int)pi(sqrt(xi));
				for(int j = i;j < B;j++){
					ret += pi(xi/primes[j]) - j;
				}
			}
			return ret;
		}
		
		public static int[] sieveEratosthenes(int n) {
			if (n <= 32) {
				int[] primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
				for (int i = 0; i < primes.length; i++) {
					if (n < primes[i]) {
						return Arrays.copyOf(primes, i);
					}
				}
				return primes;
			}

			int u = n + 32;
			double lu = Math.log(u);
			int[] ret = new int[(int) (u / lu + u / lu / lu * 1.5)];
			ret[0] = 2;
			int pos = 1;

			int[] isnp = new int[(n + 1) / 32 / 2 + 1];
			int sup = (n + 1) / 32 / 2 + 1;

			int[] tprimes = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
			for (int tp : tprimes) {
				ret[pos++] = tp;
				int[] ptn = new int[tp];
				for (int i = (tp - 3) / 2; i < tp << 5; i += tp)
					ptn[i >> 5] |= 1 << (i & 31);
				for (int j = 0; j < sup; j += tp) {
					for (int i = 0; i < tp && i + j < sup; i++) {
						isnp[j + i] |= ptn[i];
					}
				}
			}

			// 3,5,7
			// 2x+3=n
			int[] magic = { 0, 1, 23, 2, 29, 24, 19, 3, 30, 27, 25, 11, 20, 8, 4, 13, 31, 22, 28, 18, 26, 10, 7, 12, 21, 17,
					9, 6, 16, 5, 15, 14 };
			int h = n / 2;
			for (int i = 0; i < sup; i++) {
				for (int j = ~isnp[i]; j != 0; j &= j - 1) {
					int pp = i << 5 | magic[(j & -j) * 0x076be629 >>> 27];
					int p = 2 * pp + 3;
					if (p > n)
						break;
					ret[pos++] = p;
					if ((long) p * p > n)
						continue;
					for (int q = (p * p - 3) / 2; q <= h; q += p)
						isnp[q >> 5] |= 1 << q;
				}
			}

			return Arrays.copyOf(ret, pos);
		}
		
		private static class LongHashCounterL {
			public long[] keys;
			public long[] allocated;
			private int scale = 1<<2;
			private int rscale = 1<<1;
			private int mask = scale-1;
			public int size = 0;
			
			public LongHashCounterL(){
				allocated = new long[scale];
				Arrays.fill(allocated, NG);
				keys = new long[scale];
			}
			
			// if value is NG, entry is removed. (e.g. 0)
			private static final int NG = 0;
			
			public boolean containsKey(long x)
			{
				int pos = h(x)&mask;
				while(allocated[pos] != NG){
					if(x == keys[pos])return true;
					pos = pos+1&mask;
				}
				return false;
			}
			
			public long get(long x)
			{
				int pos = h(x)&mask;
				while(allocated[pos] != NG){
					if(x == keys[pos])return allocated[pos];
					pos = pos+1&mask;
				}
				return NG;
			}
			
			public long put(long x, long v)
			{
				int pos = h(x)&mask;
				while(allocated[pos] != NG){
					if(x == keys[pos]){
						long oldval = allocated[pos];
						allocated[pos] = v;
						return oldval;
					}
					pos = pos+1&mask;
				}
				if(size == rscale){
					resizeAndPut(x, v);
				}else{
					keys[pos] = x;
					allocated[pos] = v;
				}
				size++;
				return 0;
			}
			
			public long inc(long x, long v)
			{
				int pos = h(x)&mask;
				while(allocated[pos] != NG){
					if(x == keys[pos]){
						allocated[pos] += v;
						return allocated[pos];
					}
					pos = pos+1&mask;
				}
				if(size == rscale){
					resizeAndPut(x, v);
				}else{
					keys[pos] = x;
					allocated[pos] = v;
				}
				size++;
				return v;
			}
			
			public boolean remove(long x)
			{
				int pos = h(x)&mask;
				while(allocated[pos] != NG){
					if(x == keys[pos]){
						size--;
						// take last and fill rmpos
						int last = pos;
						pos = pos+1&mask;
						while(allocated[pos] != NG){
							int lh = h(keys[pos])&mask;
							// lh <= last < pos
							if(
									lh <= last && last < pos ||
									pos < lh && lh <= last ||
									last < pos && pos < lh
									){
								keys[last] = keys[pos];
								allocated[last] = allocated[pos];
								last = pos;
							}
							pos = pos+1&mask;
						}
						keys[last] = 0;
						allocated[last] = NG;
						
						return true;
					}
					pos = pos+1&mask;
				}
				return false;
			}
			
			private void resizeAndPut(long x, long v)
			{
				int nscale = scale<<1;
				int nrscale = rscale<<1;
				int nmask = nscale-1;
				long[] nallocated = new long[nscale];
				Arrays.fill(nallocated, NG);
				long[] nkeys = new long[nscale];
				for(int i = next(0);i < scale;i = next(i+1)){
					long y = keys[i];
					int pos = h(y)&nmask;
					while(nallocated[pos] != NG)pos = pos+1&nmask;
					nkeys[pos] = y;
					nallocated[pos] = allocated[i];
				}
				{
					int pos = h(x)&nmask;
					while(nallocated[pos] != NG)pos = pos+1&nmask;
					nkeys[pos] = x;
					nallocated[pos] = v;
				}
				allocated = nallocated;
				keys = nkeys;
				scale = nscale;
				rscale = nrscale;
				mask = nmask;
			}
			
			public int next(int itr)
			{
				while(itr < scale && allocated[itr] == NG)itr++;
				return itr;
			}
			
			private int h(long x)
			{
				x ^= x>>>33;
				x *= 0xff51afd7ed558ccdL;
				x ^= x>>>33;
				x *= 0xc4ceb9fe1a85ec53L;
				x ^= x>>>33;
				return (int)x;
			}
			
			@Override
			public String toString()
			{
				StringBuilder sb = new StringBuilder();
				for(int i = next(0);i < scale;i = next(i+1)){
					sb.append(",");
					sb.append(keys[i] + ":" + allocated[i]);
				}
				return sb.length() == 0 ? "" : sb.substring(1);
			}
		}
	}
	
	void run() throws Exception
	{
		is = oj ? System.in : new ByteArrayInputStream(INPUT.getBytes());
		out = new PrintWriter(System.out);
		
		long s = System.currentTimeMillis();
		solve();
		out.flush();
		tr(System.currentTimeMillis()-s+"ms");
	}
	
	public static void main(String[] args) throws Exception { new F3().run(); }
	
	private byte[] inbuf = new byte[1024];
	public int lenbuf = 0, ptrbuf = 0;
	
	private int readByte()
	{
		if(lenbuf == -1)throw new InputMismatchException();
		if(ptrbuf >= lenbuf){
			ptrbuf = 0;
			try { lenbuf = is.read(inbuf); } catch (IOException e) { throw new InputMismatchException(); }
			if(lenbuf <= 0)return -1;
		}
		return inbuf[ptrbuf++];
	}
	
	private boolean isSpaceChar(int c) { return !(c >= 33 && c <= 126); }
	private int skip() { int b; while((b = readByte()) != -1 && isSpaceChar(b)); return b; }
	
	private double nd() { return Double.parseDouble(ns()); }
	private char nc() { return (char)skip(); }
	
	private String ns()
	{
		int b = skip();
		StringBuilder sb = new StringBuilder();
		while(!(isSpaceChar(b))){ // when nextLine, (isSpaceChar(b) && b != ' ')
			sb.appendCodePoint(b);
			b = readByte();
		}
		return sb.toString();
	}
	
	private char[] ns(int n)
	{
		char[] buf = new char[n];
		int b = skip(), p = 0;
		while(p < n && !(isSpaceChar(b))){
			buf[p++] = (char)b;
			b = readByte();
		}
		return n == p ? buf : Arrays.copyOf(buf, p);
	}
	
	private char[][] nm(int n, int m)
	{
		char[][] map = new char[n][];
		for(int i = 0;i < n;i++)map[i] = ns(m);
		return map;
	}
	
	private int[] na(int n)
	{
		int[] a = new int[n];
		for(int i = 0;i < n;i++)a[i] = ni();
		return a;
	}
	
	private int ni()
	{
		int num = 0, b;
		boolean minus = false;
		while((b = readByte()) != -1 && !((b >= '0' && b <= '9') || b == '-'));
		if(b == '-'){
			minus = true;
			b = readByte();
		}
		
		while(true){
			if(b >= '0' && b <= '9'){
				num = num * 10 + (b - '0');
			}else{
				return minus ? -num : num;
			}
			b = readByte();
		}
	}
	
	private long nl()
	{
		long num = 0;
		int b;
		boolean minus = false;
		while((b = readByte()) != -1 && !((b >= '0' && b <= '9') || b == '-'));
		if(b == '-'){
			minus = true;
			b = readByte();
		}
		
		while(true){
			if(b >= '0' && b <= '9'){
				num = num * 10 + (b - '0');
			}else{
				return minus ? -num : num;
			}
			b = readByte();
		}
	}
	
	private boolean oj = System.getProperty("ONLINE_JUDGE") != null;
	private void tr(Object... o) { if(!oj)System.out.println(Arrays.deepToString(o)); }
}
