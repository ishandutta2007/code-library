

| Assumption       | Formula                     | LHS Time  | RHS Time  | LHS Space  | RHS Space|Remarks|
|:-------------    |:---------------------------:| ---------:| ---------:|-----------:|---------:|------:|
| --|--|--|--|--|--|--|
| --|$\sum_{k\le n} \gcd(k,n) = \sum_{d\|n} d\cdot\phi(n/d)$|--|--|--|--|--|
| --|Number of  squarfree nos = $\sum_{i\le \sqrt x} \mu(i)f(x/i^2)$[A057627](https://oeis.org/A057627)|--|--|--|--|--|
| $f(x)=\sum_{k\le x} \omega(k)$ [A013939](https://oeis.org/A013939) |$\sum_{i\le x} f(x/i^2)=\sum_{i\le n}\sum_{d\|i}[d \le \sqrt i]=$</br>$\sum_{i\le n}(\Omega(i) \mod 2)\cdot \lfloor{n/i}\rfloor$|$O(n^{3/4})$|--|--|--|--|
| $f(x)=\sum_{k\le x} \mu(k)$[A002321](https://oeis.org/A002321)|$\sum_{i\le x} f(x/i)=1$|$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x)=\sum_{k\le x} [\mu(k)=1]$[A070548](https://oeis.org/A070548)|$\sum_{i\le x} f(x/i) = 1 + \sum_{k \le n} \omega(k)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} [\mu(k)=-1]$[A070549](https://oeis.org/A070549)|$\sum_{i\le x} f(x/i) = \sum_{k \le n} \omega(k)$[A013939](https://oeis.org/A013939)|--|--|--|--|--|
| $f(x)=\sum_{k\le x} [\Omega(k)=1]$|$\sum_{i\le x} f(x/i) = 1 + \sum_{k \le n-1} \Omega(k)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k \cdot \mu(k)$[A068340](https://oeis.org/A068340)|$\sum_{i\le x} i \cdot f(x/i)=1$|$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x)=\sum_{k\le x} k^2 \cdot \mu(k)$|$\sum_{i\le x} i^2 \cdot f(x/i)=1$|$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x)=\sum_{k\le x} k^t \cdot \mu(k)$|$\sum_{i\le x} i^t \cdot f(x/i)=1$|$O(n^{3/4})$|$O(1)$|--|--|holds for any arbitary integer $t$|
| $f(x)=\sum_{k\le x} k \cdot \mu(k)$|$\sum_{i\le x} i^2 \cdot f(x/i)=\sum_{k\le x} k\cdot\varphi(k)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k^2 \cdot \mu(k)$|$\sum_{i\le x} i^3 \cdot f(x/i)=\sum_{k\le x} k^2\cdot\varphi(k)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k^3 \cdot \mu(k)$|$\sum_{i\le x} i^4 \cdot f(x/i)=\sum_{k\le x} k^3\cdot\varphi(k)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k^t \cdot \mu(k)$|$\sum_{i\le x} i^{t+1} \cdot f(x/i)=\sum_{k\le x} k^t\cdot\varphi(k)$|--|--|--|--|holds for any arbitary integer $t$|
| $f(x)=\sum_{k\le x} \varphi(k)$|$\sum_{i \le x}f(x/i)=x(x+1)/2$ |$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x)=\sum_{k\le x} \varphi(k)$|$\sum_{i \le x} i \cdot f(x/i)=\sum_{i\le n}\sum_{j\le n} gcd(i,j)$|$O(n^{3/4})$|--|--|--|[PE 625][GCDEX2][PE351]|
| $f(x)=\sum_{k\le x} k\cdot\varphi(k)$ [A011755](https://oeis.org/A011755)|$\sum_{i \le x} i \cdot f(x/i)=x(x+1)(2x+1)/6$|$O(n^{3/4})$|$O(1)$|--|--|[PE 448]|
| $f(x)=\sum_{k\le x} k^2\cdot\varphi(k)$|$\sum_{i \le x} i^2 \cdot f(x/i)=(x(x+1)/2)^2$|$O(n^{3/4})$|$O(1)$|--|--|--|
| --|$\sum_{m\le k\le n} lcm(k,n) =$ <br /> $(n/2) \sum_{d\|n, d > 1} d \cdot phi(d)$ <br /> $+ n - n \sum_{1\le i \le{m-1}} \frac{i}{gcd(i,n)}$|--|--|--|--|[ADDLCM]|
| f(n,k)=sum_of_coprimes_of_n_till_k|$\sum_{1\le i \le{m-1}} \frac{i}{gcd(i,n)} =$ <br /> $\sum_{d\|n, d > 1} f(\lfloor {n/d} \rfloor, \lfloor {(m-1)/d} \rfloor)$|--|--|--|--|[ADDLCM]|
| $f(x)=\sum_{k\le x} k\cdot\varphi(k)$|$\sum_{i \le x} i^2 \cdot f(x/i^2)=(x^3+3x^2+2x+3)/3$|$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x)=\sum_{k\le x} k\cdot\varphi(k)$|$\sum_{i \le x} i^3 \cdot f(x/i^3)=f(x)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k\cdot\varphi(k)$|$\sum_{i \le x} i^4 \cdot f(x/i^4)=f(x)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k\cdot\varphi(k)$|$\sum_{i \le x} i^t \cdot f(x/i^t)=f(x)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k^2\cdot\varphi(k)$|$\sum_{i \le x} i^3 \cdot f(x/i^3)=f(x)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k^2\cdot\varphi(k)$|$\sum_{i \le x} i^4 \cdot f(x/i^4)=f(x)$|--|--|--|--|--|
| $f(x)=\sum_{k\le x} k^2\cdot\varphi(k)$|$\sum_{i \le x} i^t \cdot f(x/i^t)=f(x)$|--|--|--|--|--|
| --|$\sum_{i\le n}\sum_{j\le n} gcd(i,j) =$ <br /> $\sum_{d\le n} \mu(d)\lfloor n/d \rfloor ^2$|--|--|--|--|[PE 625][GCDEX2]|
| --|$\sum_{k\le n} k\cdot \varphi(k) \ \ =$ <br /> $\sum_{d\le n}{\mu(d)\cdot d \cdot S\left(\left[\frac{n}{d}\right]\right)}, \tag{1}$  <br /> where $S(i)=\sum i^2$|--|--|--|--|--|
| $f(x)=\sum_{i\le x} \sum_{j\le x} \sigma_2(i\cdot j)$|$\sum_{i \le x} i^2 \cdot f(x/i) = (\sum_{k \le n} k^2 \cdot \lfloor {n/k} \rfloor)^2 $|--|--|--|--|--|
| $f(x)=\sum_{i\le x} \sum_{j\le x} \sigma_1(i\cdot j)$|$\sum_{i \le x} i \cdot f(x/i) = (\sum_{k\le n} k\cdot \lfloor {n/k} \rfloor)^2 $|$O(n^{3/4})$|$O(n^{1/2})$|--|--|[PE 439]|
| $f(x)=\sum_{i\le x} \sum_{j\le x} \sigma_0(i\cdot j)$|$\sum_{i \le x} f(x/i) = (\sum_{k\le n} \lfloor {n/k} \rfloor)^2 $|--|--|--|--|--|--|
| $\phi(x, a)=\sum_{i\le x} \prod_{j\le a} \[i \mod prime_{j} \ne 0]$<br/>$\phi(x,a)=$ Legendre's prime counting formula|$\sum_{i \le x}\phi(\lfloor {x/prime_{i}} \rfloor, i-1) = i$|$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x, a) = \sum_{i\le x} \prod_{j\le a} [i \mod prime_{j} \ne 0]$|$\sum_{i \le x}  prime_{i} \cdot f(\lfloor {x/prime_{i}} \rfloor, i-1) = \sum_{i \le x} spf(i)$<br/>where spf = smallest prime factor|$O(n^{3/4})$|--|--|--|--|
| $f(x, a) = \sum_{i\le x} i \cdot \prod_{j\le a}[i \mod prime_{j} \ne 0]$|$\sum_{i \le x}  prime_{i} \cdot f(\lfloor {x/prime_{i}} \rfloor, i-1) = x(x+3)/2$|$O(n^{3/4})$|$O(1)$|--|--|--|
| $f(x, a)= \sum_{i\le x} i \cdot \prod_{j\le a} [i \mod prime_{j} \ne 0]$|$\sum_{i \le x} f(\lfloor {x/prime_{i}} \rfloor, i-1) = \sum_{2 \le k \le x} k/spf(k)$|--|--|--|--|--|
| $f(x)=\sum_{i\le x} \sigma_0(i)$[A006218](https://oeis.org/A006218)|$\sum_{i \le x} f(x/i) = \(tau<=)_3(n)$[A061201](https://oeis.org/A061201)|$O(n^{3/4})$|$O(n^{1/4})$|--|--|--|
| --|$\sum_{i\le x} \sigma_0(i) = \sum_{k \le n} \lfloor{n/k}\rfloor =$</br>$\sum_{k \le \sqrt n} \lfloor{n/k}\rfloor + \sum_{v \le \sqrt n} v\cdot(\lfloor{n/v}\rfloor-\lfloor{n/(v+1)}\rfloor)$|--|$O(n^{1/2})$|--|--|--|
| --|$\sum_{i\le x} \sigma_0(i) = hyperbolic\ method$|--|$O(n^{1/3})$|--|--|--|
| $f(x)=\sum_{i\le x} \sigma_0(i^2)$[A061503](https://oeis.org/A061503)|$\sum_{i \le x} f(x/i) = \sum_{k \le n} \sigma_0(k)^2$[A061502](https://oeis.org/A061502)|$O(n^{3/4})$|--|--|--|--|
| --|$\sum_{i\le x} \sigma_0(i)=$<br />$\small{\sum_{1 \le k \le \lfloor \sqrt n \rfloor} (2.\sum_{1 \le j \le \lfloor(\sqrt{n/k^2})\rfloor} \lfloor{n/(j.k^2)}\rfloor - \lfloor(\sqrt{n/k^2})^2\rfloor)}$|--|--|--|--|--|
| $f(x)=\sum_{i\le x} \sigma_1(i)$[A024916](https://oeis.org/A024916)|$\sum_{i\le x} f(x/i)= \sum_{d\|x} \sigma_1(d)$[A280077](https://oeis.org/A280077)</br>$= \sum_{k\le x}\sigma_1(k)\lfloor{x/k}\rfloor$|$O(n^{3/4})$|$O(n^{1/4})$|--|--|--|
| $f(x)=\sum_{i\le x} \sigma_2(i)$[A064602](https://oeis.org/A064602)|$\sum_{i\le x} f(x/i)= \sum_{k\le x}\sigma_2(k)\lfloor{x/k}\rfloor$|--|--|--|--|--|
| --|--|--|--|--|--|--|

-  
TODO
152
243
319
338
347

360
370
379
388
415
432
465
556

