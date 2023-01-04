#include <bits/stdc++.h>
using namespace std;
//  This code uses min heap instead of solving pells equation
//  This code is 100x faster than pythons(for p=47)
//  https://en.wikipedia.org/wiki/St%C3%B8rmer%27s_theorem
std::vector<std::pair<long long, long long>>
stormer(std::vector<unsigned int> primes) {
  std::reverse(primes.begin(), primes.end());
  std::vector<std::pair<long long, long long>> solutions;
  // remove large prime depending on user input
  while (!primes.empty() && primes.front() > primes[0])
    primes.erase(primes.begin());

  // min-heap
  std::priority_queue<unsigned long long, std::vector<unsigned long long>,
                      std::greater<unsigned long long>> next;
  // "seed" value
  next.push(1);

  // will store the final result
  unsigned long long sum = 0;

  // 47-smooth number from previous iteration
  unsigned long long last = 1; // must not be zero, any other value is fine

  // I saw in my first runs that no solutions were found beyond the 10 millionth
  // 47-smooth number
  for (unsigned int iteration = 0; iteration < 10000000; iteration++) {
    // fetch next 47-smooth number
    auto current = next.top();
    next.pop();

    // two consecutive 47-smooth numbers ? => T(last) is 47-smooth, too
    if (last + 1 == current) {
      solutions.push_back(std::pair<long long, long long>(last, last + 1));
      sum += last;
      // std::cout << "T(" << last << "), sum=" << sum << " @ iteration " <<
      // iterations << " " << next.size() << std::endl;
    }

    // remember for next iteration
    last = current;

    // find further 47-smooth numbers
    for (auto p : primes) {
      auto todo = current * p;

      // heuristic: ignore numbers beyond the largest 47-smooth T(n)
      if (todo < 1111111111111ULL)
        next.push(todo);

      // any prime smaller than the largest prime yields a number already found
      // in "next"
      if (current % p == 0)
        break;
    }
  }

  // std::cout << sum << std::endl;
  return solutions;
}

int main() {
  std::vector<unsigned int> primes = {2,  3,  5,  7,  11, 13, 17,
                                      19, 23, 31, 37, 41, 43, 47};
  std::vector<std::pair<long long, long long>> solutions = stormer(primes);
  for (auto solution : solutions)
    cout << solution.first << " " << solution.second << endl;
  return 0;
}
