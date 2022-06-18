#include <bits/stdc++.h>
using namespace std;

//длина, до которой эффективнее обычное умножение
#define EFF_SIZE 1024

//обратные к n = 2^k по модулю 2^n + 1 для k = 1, 2, 3, 4, 5, 6
// long invn[7] = {0, 3, 13, 225, 61441, 4160749569, 18158513697557839873};
//наименьшая степень двойки, превосходящая a
int exp2ab(int a) {
  int k = 1;
  while (a > k) {
    k <<= 1;
  }
  return k;
}
int log2(int a) {
  int k = 0;
  while (a > 1) {
    k++;
    a >>= 1;
  }
  return k;
}

// 2^n
int exp2(int n) {
  if (!n)
    return 1;
  else
    return 1 << n;
}

unsigned int bit_inverse(unsigned int a, int n) {
  int Shift = n - 1;
  unsigned int LowMask = 1;
  unsigned int HighMask = 1 << Shift;
  unsigned int R;
  for (R = 0; Shift >= 0; LowMask <<= 1, HighMask >>= 1, Shift -= 2)
    R |= ((a & LowMask) << Shift) | ((a & HighMask) >> Shift);
  return R;
}

//перестановка входного массива
template <class T> vector<T> trans(vector<T> x) {
  int n = log2(x.size());
  vector<T> tmp;
  tmp = x;
  for (int i = 0; i < exp2(n); i++) {
    x[i] = tmp[bit_inverse(i, n)];
  }
  return x;
}

class longint {
public:
  vector<char> number;
  longint operator=(const longint &inumber) {
    number = inumber.number;
    return *this;
  }
  bool operator>(const longint &inumber) {
    int i;
    if (number.size() > inumber.number.size())
      return true;
    if (number.size() < inumber.number.size())
      return false;
    for (i = number.size() - 1; i >= 0; i--) {
      if (number[i] > inumber.number[i])
        return true;
      if (number[i] < inumber.number[i])
        return false;
    }
    return false;
  }
  bool operator==(const longint &inumber) { return number == inumber.number; }
  longint() { number = vector<char>(1, 0); }
  longint(char *filename) {
    FILE *in = fopen(filename, "r");
    char buf;
    int read = fscanf(in, "%c", &buf);
    while (read != -1) {
      number.push_back(buf - 48);
      read = fscanf(in, "%c", &buf);
    }
    fclose(in);
  }
  longint(const longint &inumber) { number = inumber.number; }
  longint(const vector<char> &ivec) { number = ivec; }
  longint(long numb) {
    if (numb == 0)
      number.push_back(0);
    while (numb > 0) {
      number.push_back(numb & 1);
      numb >>= 1;
    }
  }
  void print() {
    for (int i = number.size() - 1; i >= 0; i--) {
      printf("%c", number[i] + 48);
    }
    printf("\n");
  }
  void print(char *str) {
    // FILE* out = fopen(str, "w");
    for (int i = number.size() - 1; i >= 0; i--) {
      printf("%c", number[i] + 48);
    }
    printf("\n");
    // fclose(out);
  }
  bool iszero() { return (number.size() == 1) && (number[0] == 0); }
  size_t size() { return number.size(); }
  //операции сдвига с сохранением числа разрядов
  void lshift() {
    number.pop_back();
    number.insert(number.begin(), 0);
  }
  void lshift(int k) {
    int i;
    for (i = 0; i < k; i++) {
      number.pop_back();
      number.insert(number.begin(), 0);
    }
  }
  void rshift() {
    number.erase(number.begin());
    number.push_back(0);
  }
  void rshift(int k) {
    int i;
    for (i = 0; i < k; i++) {
      number.erase(number.begin());
      number.push_back(0);
    }
  }
  //операции умножения и целочисленного деления на два
  void mult2() { number.insert(number.begin(), 0); }
  void mult2(int k) { number.insert(number.begin(), k, 0); }
  void div2() {
    number.erase(number.begin());
    if (number.empty())
      number = vector<char>(1, 0);
  }
  void div2(int k) {
    int m = min((int)number.size(), k);
    number.erase(number.begin(), number.begin() + m);
    if (number.empty())
      number = vector<char>(1, 0);
  }
  //сложение
  void add(const longint inumber) {
    if (inumber.number.size() == 1 && inumber.number[0] == 0)
      return;
    vector<char> result(max(number.size(), inumber.number.size()) + 1);
    vector<char> minlength =
        (number.size() > inumber.number.size()) ? inumber.number : number;
    vector<char> maxlength =
        (number.size() <= inumber.number.size()) ? inumber.number : number;
    int min_length = (number.size() > inumber.number.size())
                         ? inumber.number.size()
                         : number.size();
    int max_length = (number.size() <= inumber.number.size())
                         ? inumber.number.size()
                         : number.size();
    int i;
    char bitcarry = 0;
    for (i = 0; i < min_length; i++) {
      result[i] = minlength[i] ^ maxlength[i] ^ bitcarry;
      bitcarry = (minlength[i] ^ maxlength[i]) & bitcarry |
                 (minlength[i] & maxlength[i]);
    }
    for (i = min_length; i < result.size() - 1; i++) {
      result[i] = maxlength[i] ^ bitcarry;
      bitcarry = maxlength[i] & bitcarry;
    }
    result[result.size() - 1] = bitcarry;
    number = result;
    while (number[number.size() - 1] == 0 && number.size() > 1) {
      number.pop_back();
    }
  }
  //добавление числа, сдвинутого на k разрядов влево
  void add(const longint inumber, int k) {
    if (inumber.number.size() == 1 && inumber.number[0] == 0)
      return;
    vector<char> result(max(number.size(), inumber.number.size()) + k + 1);

    vector<char> first = number;
    vector<char> second = inumber.number;

    if (first.size() < result.size())
      first.insert(first.end(), result.size() - first.size(), 0);
    if (second.size() < result.size())
      second.insert(second.end(), result.size() - second.size() - k, 0);

    int i;
    char bitcarry = 0;
    for (i = 0; i < k; i++)
      result[i] = first[i];
    for (i = k; i < result.size(); i++) {
      result[i] = first[i] ^ second[i - k] ^ bitcarry;
      bitcarry =
          (first[i] ^ second[i - k]) & bitcarry | (first[i] & second[i - k]);
    }
    number = result;
    while (number[number.size() - 1] == 0 && number.size() > 1) {
      number.pop_back();
    }
  }
  //вычитание
  void sub(const longint inumber) {

    vector<char> result(max(number.size(), inumber.number.size()));
    vector<char> b = inumber.number;
    vector<char> a = number;
    int b_length = inumber.number.size();
    int a_length = number.size();
    int i;
    char bitcarry = 0;
    for (i = 0; i < b_length; i++) {
      result[i] = b[i] ^ a[i] ^ bitcarry;
      bitcarry =
          (bitcarry ^ 1) & (a[i] ^ 1) & b[i] | bitcarry & (a[i] ^ 1 | b[i]);
    }
    for (i = b_length; i < result.size(); i++) {
      result[i] = a[i] ^ bitcarry;
      bitcarry = (a[i] ^ 1) & bitcarry;
    }
    number = result;
    while (number[number.size() - 1] == 0 && number.size() > 1) {
      number.pop_back();
    }
  }
  //противоположный элемент по модулю 2^n+1, где n - количество бит в записи
  //числа
  void inv(int n) {
    if (number.size() == 1 && number[0] == 0)
      return;
    longint module(1);
    module.mult2(n);
    module.add(longint(1));

    longint temp(number);

    while (temp > module)
      temp.sub(module);

    module.sub(temp);
    number = module.number;
  }

  longint operator+(const longint &inumber) {
    longint result(number);
    result.add(inumber);
    return result;
  }
  longint operator-(const longint &inumber) {
    longint result(number);
    result.sub(inumber);
    return result;
  }

  //отбросить старшие разряды, начиная с n
  void lowbits(int n) {
    if (number.size() > n)
      number.erase(number.begin() + n, number.end());
  }
  //удаление незначащих нулей
  void delete_zeros() {
    while (number[number.size() - 1] == 0 && (number.size() > 1))
      number.pop_back();
  }
};
unsigned int toint(longint a) {
  if (a.size() > 32)
    return -1;
  int result = 0;
  int i;
  for (i = 0; i < a.size(); i++) {
    result += (1 << i) * a.number[i];
  }
  return result;
}
//обычное умножение
longint mult(longint a, longint b) {
  int i;
  longint result;
  longint temp = b;
  for (i = 0; i < a.number.size(); i++) {
    if (a.number[i] == 1) {
      result.add(temp, i);
    }
  }
  return result;
}
longint mod(longint inumber, int n, int r) {

  // temp = 2^n + 1
  inumber.delete_zeros();
  vector<char> temp(n + 1, 0);
  temp[n] = 1;
  temp[0] = 1;
  longint module(temp);

  r %= 2 * n;
  if (r == 0) {
    if (module > inumber)
      return inumber;
    if (inumber == module)
      return longint((long)0);
  }

  inumber.mult2(r);
  if (module > inumber)
    return inumber;

  int i;
  longint result;
  longint sum1((long)0);
  longint sum2((long)0);

  while (inumber.size() % n != 0)
    inumber.number.push_back(0);
  for (i = 0; i < inumber.size() / n; i++) {
    if (i % 2 == 0) {
      sum1.add(longint(vector<char>(inumber.number.begin() + i * n,
                                    inumber.number.begin() + (i + 1) * n)));
    } else {
      sum2.add(longint(vector<char>(inumber.number.begin() + i * n,
                                    inumber.number.begin() + (i + 1) * n)));
    }
  }
  while (sum2 > sum1)
    sum1.add(module);

  inumber = sum1 - sum2;
  inumber.delete_zeros();
  return mod(inumber, n, 0);
}
longint mod(longint inumber, int n) { return mod(inumber, n, 0); }
//ДПФ в кольце Z_2^n+1
vector<longint> rdft(vector<longint> &a, int n) {

  vector<longint> result(a.size(), longint((long)0));
  int K = a.size();
  // 2^n + 1
  longint temp;
  int i, j;
  int l, q, r;
  int N = a.size();
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      l = 2 * n * i * j / K;
      q = l / n;
      r = l % n;
      temp = a[j];
      temp.mult2(r);
      temp = mod(temp, n);
      if (q & 1) {
        temp.inv(n);
      }
      result[i].add(temp);
      result[i] = mod(result[i], n);
    }
  }
  return result;
}

//ОДПФ в кольце Z_2^n+1
vector<longint> irdft(vector<longint> &a, int n) {

  vector<longint> result(a.size(), longint());
  int K = a.size();
  longint temp;
  int i, j;
  int l, q, r;
  int N = a.size();
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      l = 2 * n * i * j - 2 * n * i * j / K;
      q = l / n;
      r = l % n;
      temp = a[j];
      temp.mult2(r);
      temp = mod(temp, n);
      if (q & 1) {
        temp.inv(n);
      }
      result[i].add(temp);
      result[i] = mod(result[i], n);
    }
  }
  //умножение на обратный к K по модулю 2^n + 1
  int k = log2(K);
  for (i = 0; i < N; i++) {
    result[i].mult2(2 * n - k);
    result[i] = mod(result[i], n);
  }
  return result;
}
//БПФ в кольце Z_2^n+1
vector<longint> rfft(vector<longint> &a, int n) {
  int K = (int)a.size();
  int N = log2(K);
  a = trans(a);

  int i, j, k;

  longint module(1);
  module.mult2(n);
  module.add(1);

  longint u, v;

  for (i = 2; i <= K; i <<= 1) {
    int root = 2 * n / i;
    for (j = 0; j < K; j += i) {
      int omega = 0;
      for (int k = 0; k < i / 2; k++) {
        u = a[j + k];
        v = a[j + k + i / 2];
        v = mod(v, n, omega);

        a[j + k] = u + v;
        a[j + k + i / 2] = module + module + u - v;
        omega += root;

        a[j + k] = mod(a[j + k], n);
        a[j + k + i / 2] = mod(a[j + k + i / 2], n);
      }
    }
  }

  return a;
}
vector<longint> irfft(vector<longint> &a, int n) {
  int K = (int)a.size();
  int N = log2(K);
  a = trans(a);

  int i, j, k;

  vector<char> temp(n + 1, 0);
  temp[n] = 1;
  temp[0] = 1;
  longint module(temp);

  longint u, v;

  for (i = 2; i <= K; i <<= 1) {
    // printf("\n%d\n", i);
    int root = 2 * n - 2 * n / i;
    for (j = 0; j < K; j += i) {
      int omega = 0;
      for (int k = 0; k < i / 2; k++) {
        u = a[j + k];
        v = a[j + k + i / 2];
        v = mod(v, n, omega);

        a[j + k] = u + v;
        a[j + k + i / 2] = module + module + u - v;
        omega += root;

        a[j + k] = mod(a[j + k], n);
        a[j + k + i / 2] = mod(a[j + k + i / 2], n);
      }
    }
  }
  //умножение на обратный к K по модулю 2^n + 1
  for (i = 0; i < K; i++) {
    v = mod(a[i], n, 2 * n - N);
    if ((2 * n - N) < 0) {
      printf("error");
    }
  }
  return a;
}
//разбиение на фрагменты длины block_size
vector<longint> decomposition(longint number, int block_size) {
  vector<longint> decomposed_vector;
  int i = 0;
  while ((i + 1) * block_size <= number.size()) {
    decomposed_vector.push_back(
        longint(vector<char>(number.number.begin() + block_size * i,
                             number.number.begin() + block_size * (i + 1))));
    i++;
  }
  return decomposed_vector;
}
//взвешивание
vector<longint> weight(vector<longint> in_vec, int n) {
  int K = in_vec.size();
  int i, l, q, r;
  for (i = 0; i < K; i++) {
    l = n * i / K;
    q = l / n;
    r = l % n;

    in_vec[i] = mod(in_vec[i], n, r);
    if (q & 1) {
      in_vec[i].inv(n);
    }
    in_vec[i] = mod(in_vec[i], n);
  }
  return in_vec;
}
//обратное взвешивание
vector<longint> unweight(vector<longint> in_vec, int n) {
  int K = in_vec.size();
  int i, l, q, r;
  for (i = 0; i < K; i++) {
    l = 2 * n * i - n * i / K;
    q = l / n;
    r = l % n;

    in_vec[i] = mod(in_vec[i], n, r);
    if (q & 1) {
      in_vec[i].inv(n);
    }
    in_vec[i] = mod(in_vec[i], n);
  }
  return in_vec;
}
//перенос
vector<longint> carry(vector<longint> in_vec, int block_size) {
  longint temp;
  vector<longint> result(in_vec);
  result.push_back(longint(vector<char>(block_size, 0)));
  int N = in_vec.size();
  int i;
  for (i = 0; i < N; i++) {
    if (result[i].size() < block_size)
      result[i].number.insert(result[i].number.end(),
                              block_size - result[i].size(), 0);
  }
  for (i = 0; i < N; i++) {
    temp = result[i];
    if (result[i].size() > block_size)
      result[i + 1].add(longint(
          vector<char>(temp.number.begin() + block_size, temp.number.end())));
    result[i].lowbits(block_size);
  }
  return result;
}
//"склейка"
longint composition(vector<longint> vec) {
  vector<char> result;
  int i;
  for (i = 0; i < vec.size(); i++) {
    result.insert(result.end(), vec[i].number.begin(), vec[i].number.end());
  }
  return longint(result);
}
//вывод
void print(vector<longint> in_vec) {
  int i;
  for (i = 0; i < in_vec.size(); i++) {
    in_vec[i].print();
  }
}
vector<longint> neg_conv(vector<longint> a, vector<longint> b, int n) {

  int N = a.size();
  vector<longint> result(N);
  a = weight(a, n);
  b = weight(b, n);

  a = rfft(a, n);
  b = rfft(b, n);
  for (int i = 0; i < N; i++) {
    result[i] = mult(a[i], b[i]);
    result[i] = mod(result[i], n);
  }
  result = irfft(result, n);
  result = unweight(result, n);
  return result;
}
//умножение методом Шёнхаге-Штрассена
longint ss_mult(longint a, longint b, int N) {
  int min_size = min(a.number.size(), b.number.size());
  int max_size = max(a.number.size(), b.number.size());
  longint result((long)0);
  if (min_size <= EFF_SIZE) {
    // printf("return usual mult, N = %d\n", N);
    return mod(mult(a, b), N);
  }
  int i, j;
  if (max_size > N) {
    a = mod(a, N);
    b = mod(b, N);
    return ss_mult(a, b, N);
  }

  int K = sqrt((double)N);
  K = exp2ab(K);
  int k = log2(K);
  int M = N / K;
  int n = 2 * M + k;
  n = exp2ab(n);
  if (n < K)
    n = K;
  if (n >= N) {
    return mod(mult(a, b), N);
  }

  //дополняем числа нулями
  while (a.size() < N)
    a.number.push_back(0);
  while (b.size() < N)
    b.number.push_back(0);

  //разбиваем на слова
  vector<longint> a_dec, b_dec;
  // printf("dec\n");
  a_dec = decomposition(a, M);
  b_dec = decomposition(b, M);

  //отрицательно обёрнутая свёртка

  vector<longint> result_dec(K);

  // printf("wei\n");
  a_dec = weight(a_dec, n);
  b_dec = weight(b_dec, n);

  // printf("fft\n");
  a_dec = rfft(a_dec, n);
  b_dec = rfft(b_dec, n);

  // printf("mul\n");
  for (int i = 0; i < K; i++) {
    result_dec[i] = ss_mult(a_dec[i], b_dec[i], n);
  }

  result_dec = irfft(result_dec, n);
  // printf("unwei\n");
  result_dec = unweight(result_dec, n);
  //перенос
  // printf("car\n");
  result_dec = carry(result_dec, M);

  //склейка
  result = composition(result_dec);
  result.delete_zeros();
  result = mod(result, N);
  return result;
}
longint ss_mult(longint a, longint b) {
  int min_size = min(a.number.size(), b.number.size());
  longint result((long)0);
  if (min_size < EFF_SIZE)
    return mult(a, b);
  else {
    int N = a.number.size() + b.number.size();
    N = exp2ab(N);
    result = ss_mult(a, b, N);
  }
  return result;
}

int main() {
  longint a;
  longint b;
  int i, j;

  longint result;
  vector<char> temp(3, 1);
  a = longint(temp);
  b = longint(temp);
  a.add(temp, 1);
  a.print();

  vector<long> time_ord_mult;
  vector<long> time_ss_mult;
  vector<long> lengths;

  long time;
  int number_of_tests = 10;

  for (i = 1024; i < 1024 + 100 * 3; i += 100) {
    lengths.push_back(i);
    temp.resize(i, 1);
    a = longint(temp);
    b = longint(temp);

    printf("Length: %d\n", i);

    // time = GetTickCount();
    for (j = 0; j < number_of_tests; j++) {
      result = mult(a, b);
    }
    result.print();
    // time = GetTickCount() - time;
    // printf("Time_ord: %lf\n", time/(number_of_tests*1000.0));
    // time_ord_mult.push_back(time/number_of_tests);

    // time = GetTickCount();
    for (j = 0; j < number_of_tests; j++) {
      result = ss_mult(a, b);
    }
    result.print();
    // time = GetTickCount() - time;
    // printf("Time_ss: %lf\n", time/(number_of_tests*1000.0));
    time_ss_mult.push_back(time / number_of_tests);
  }

  // // FILE *out = fopen("time_ord", "w+");
  // for(i = 0; i < time_ord_mult.size(); i++)
  // 	printf("%d ", time_ord_mult[i]);
  // // fclose(out);

  // // out = fopen("time_ss", "w+");
  // for(i = 0; i < time_ss_mult.size(); i++)
  // 	printf("%d ", time_ss_mult[i]);
  // // fclose(out);

  // out = fopen("lengths", "w+");
  for (i = 0; i < lengths.size(); i++)
    // printf("%d ", lengths[i]);
    cout << lengths[i] << " ";
  // fclose(out);

  getchar();

  return 0;
}
