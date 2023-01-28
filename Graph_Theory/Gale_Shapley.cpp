// An example about Gale–Shapley algorithm
//
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <numeric>

using namespace std;
using ll = long long;

namespace IO {
    template<typename T>
    inline
    void read(T &t) {
        int n = 0;
        int c = getchar();
        t = 0;
        while (!isdigit(c)) n |= c == '-', c = getchar();
        while (isdigit(c)) t = t * 10 + c - 48, c = getchar();
        if (n) t = -t;
    }

    template<typename T, typename... Args>
    inline
    void read(T &t, Args &... args) {
        read(t);
        read(args...);
    }

    template<typename T>
    inline void write(T x) {
        if (x < 0) x = -x, putchar('-');
        if (x > 9) write(x / 10);
        putchar(x % 10 + 48);
    }

    template<typename T>
    inline void writeln(T x) {
        write(x);
        putchar('\n');
    }
}

namespace StableMatching {
    vector<queue<int> > m_pref;
    vector<vector<int> > w_pref;
    vector<int> engaged;    // w to m
    queue<int> free_man;
    vector<bool> is_woman_free;

    void reset() {
        vector<queue<int> >().swap(m_pref);
        vector<vector<int> >().swap(w_pref);
        vector<int>().swap(engaged);
        queue<int>().swap(free_man);
        vector<bool>().swap(is_woman_free);
    }

    void init(int n) {
        m_pref.resize(n + 1);   // start from 1
        w_pref.resize(n + 1, vector<int>(n + 1));
        engaged.resize(n + 1);
        for (int i = 1; i <= n; ++i) {
            free_man.push(i);
        }
        is_woman_free.resize(n + 1, true);
    }

    void GaleShapley() {
        while (!free_man.empty()) {
            int m = free_man.front();
            if (m_pref.at(m).empty()) {
                free_man.pop();
                continue;
            }
            int w = m_pref.at(
                    m).front();   // first woman on m's list to whom m has not yet proposed
            if (is_woman_free.at(w)) {
                engaged.at(w) = m;
                m_pref.at(m).pop();
                free_man.pop();
                is_woman_free.at(w) = false;
            } else {
                int m_ = engaged.at(w);
                if (w_pref.at(w).at(m) >
                    w_pref.at(w).at(m_)) { // w prefers m to m'
                    free_man.push(m_);
                    engaged.at(w) = m;
                    free_man.pop();
                }
                m_pref.at(m).pop();
            }
        }
    }
}

int main() {

    int t;
    IO::read(t);
    while (t--) {
        int n, id, temp;
        IO::read(n);
        StableMatching::init(n);
        // women
        for (int j = 0; j < n; ++j) {
            IO::read(id);
            for (int preference = n, i = 0; i < n; --preference, ++i) {
                IO::read(temp);
                StableMatching::w_pref.at(id).at(temp) = preference;
            }
        }
        // men
        for (int j = 0; j < n; ++j) {
            IO::read(id);
            for (int i = 0; i < n; ++i) {
                IO::read(temp);
                StableMatching::m_pref.at(id).push(temp);
            }
        }
        StableMatching::GaleShapley();
        // sort and output
        vector<pair<int, int> > output(
                StableMatching::engaged.size()); // m to w
        for (int j = 1; j <= n; ++j) {
            output.at(j) = make_pair(StableMatching::engaged.at(j), j);
        }
        sort(output.begin(), output.end());
        for (int j = 1; j <= n; ++j) {
            IO::write(output.at(j).first);
            putchar(' ');
            IO::writeln(output.at(j).second);
        }
        StableMatching::reset();
    }
}
//Luogu/SP243