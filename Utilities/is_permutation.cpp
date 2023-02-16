// https://stackoverflow.com/a/3219151/20438572
inline bool is_permutation_fast(ull a, ull b){
    int c[10]={0,0,0,0,0,0,0,0,0,0};
    while (a) { c[a%10]++; a/=10; }
    while (b) { c[b%10]--; b/=10; }
    int res=1;
    for (int i=0;i<10;i++) res &= c[i]==0;
    return(res?1:0);
}