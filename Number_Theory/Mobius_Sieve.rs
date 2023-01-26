//calcules till 10^7 in 0.1 sec
//1000x faster than c++ version

fn compute_mu(max: usize) -> Vec<isize> {
    let mut mu = vec![0; max + 1];
    let mut primes = vec![];
    let mut sieve = vec![false; max + 100];
    mu[1] = 1;
    for i in 2..=max {
        if !sieve[i as usize] {
            primes.push(i);
            mu[i] = -1;
        }
        for p in &primes {
            let p = *p;
            if i * p > max {
                break;
            }
            sieve[i * p] = true;
            if i % p == 0 {
                mu[i * p] = 0;
                break;
            } else {
                mu[i * p] = -mu[i];
            }
        }
    }

    mu
}

fn main() {
    let n = 100_000_000_000_000usize;
    let nsq = (n as f64).sqrt() as usize;
    let mu = compute_mu(nsq);
    let mut result = 0isize;
    for s in 1..=10 {
        println!("{}", mu[s]);
    }
}


