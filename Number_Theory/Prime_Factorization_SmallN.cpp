vector<int> factorize(int n) {
    vector<int> res;
    for (int i = 2; n > 1 && i < 1000; i += 1 + (i % 2))
        if (n % i == 0) {
            res.push_back(i);
            while(n % i == 0)
                n /= i;
        }
    if (n > 1)
        res.push_back(n);
    return res;
}
