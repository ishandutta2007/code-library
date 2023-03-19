
// Supports till 10^7 (to be precise till 2*3*5*7*11*13*17*19)
ll sum_of_coprimes_of_n_till_k(int n, int k) {
    auto primes = get_unique_sorted_prime_factors(n);
    int ans = T_mod(k);
    ll pp[7];
    for(int i0 = 0; i0 < primes.size(); i0++) {
        ll p1 = primes[i0];
        pp[1] = p1;
        if (pp[1] > k) break;
        // # we have to remove the series of p under k ie p+2p+3p....[k/p]p
        ll t1 = k / pp[1];
        // cout << "pp1=" << pp1 << "k=" << k << "terms=" << t1 << endl;
        ans = (ans - mulmod(pp[1], T_mod(t1), MOD) + MOD) % MOD;
        for(int i1 = i0 + 1; i1 < primes.size(); i1++) {
            ll p2 = primes[i1];
            pp[2] = pp[1] * p2;
            if (pp[2] > k) break;
            // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
            ll t2 = k / pp[2];
            // cout << "pp2=" << pp2 << "k=" << k << "terms=" << t2 << endl;
            ans = (ans + mulmod(pp[2], T_mod(t2), MOD) ) % MOD;
            for(int i2 = i1 + 1; i2 < primes.size(); i2++) {
                ll p3 = primes[i2];
                pp[3] = pp[2] * p3;
                if (pp[3] > k) break;
                // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                ll t3 = k / pp[3];
                // cout << "pp3=" << pp3 << "k=" << k << "terms=" << t3 << endl;
                ans = (ans - mulmod(pp[3], T_mod(t3), MOD)+ MOD) % MOD;
                for(int i3 = i2 + 1; i3 < primes.size(); i3++) {
                    ll p4 = primes[i3];
                    pp[4] = pp[3] * p4;
                    if (pp[4] > k) break;
                    // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                    ll t4 = k / pp[4];
                    // cout << "pp4=" << pp4 << "k=" << k << "terms=" << t4 << endl;
                    ans = (ans + mulmod(pp[4], T_mod(t4), MOD) ) % MOD;
                    for(int i4 = i3 + 1; i4 < primes.size(); i4++) {
                        ll p5 = primes[i4];
                        pp[5] = pp[4] * p5;
                        if (pp[5] > k) break;
                        // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                        ll t5 = k / pp[5];
                        // cout << "pp5=" << pp5 << "k=" << k << "terms=" << t5 << endl;
                        ans = (ans - mulmod(pp[5],T_mod(t5), MOD)+ MOD) % MOD;
                        // for (int i5=i4 + 1; i5 < primes.size(); i5++) {
                        //     ll p6 = primes[i5];
                        //     ll pp6 = pp5 * p6;
                        //     if (pp6 > k) break;
                        //     // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                        //     ll t6 = k / pp6;
                        //     // cout << "pp6=" << pp6 << "k=" << k << "terms=" << t6 << endl;
                        //     ans = (ans + mulmod(pp6, T_mod(t6), MOD)) % MOD;
                        //     // for(int i6 = i5 + 1; i6 < primes.size(); i6++) {
                        //     //     ll p7 = primes[i6];
                        //     //     ll pp7 = pp6 * p7;
                        //     //     if (pp7 > k) break;
                        //     //     // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                        //     //     ll t7 = k / pp7;
                        //     //     // cout << "pp7=" << pp7 << "k=" << k << "terms=" << t7 << endl;
                        //     //     ans = (ans - mulmod(pp7, T_mod(t7), MOD) + MOD) % MOD;
                        //     //     // for (int i7=i6 + 1; i7 < primes.size(); i7++) {
                        //     //     //     ll p8 = primes[i7];
                        //     //     //     ll pp8 = pp7 * p8;
                        //     //     //     if (pp8 > k) break;
                        //     //     //     // # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                        //     //     //     ll t8 = k / pp8;
                        //     //     //     // cout << "pp8=" << pp8 << "k=" << k << "terms=" << t8 << endl;
                        //     //     //     ans = (ans + pp8 * T_mod(t8)) % MOD;
                        //     //     // }
                        //     // }
                        // }
                    }
                }
            }
        }
    }
    return (ll)ans;
}
