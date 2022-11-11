// I am trying to iterate in ascending order through numbers
// that only consist of digits 2, 3, 5 and 7. However the
// initial input number may contain other digits.
// But after the first iteration we will be dealing
// with strictly the digits 2, 3, 5 and 7 only.

int generate_next_number(int num, static int mapdigit[] = {2, 2, 3, 5, 5, 7, 7,
                                                           2, 2, 2}) {
  int digit = num % 10;
  int nextdigit = mapdigit[digit];
  num /= 10;
  if (nextdigit < digit)
    num = generate_next_number(num);
  return num * 10 + nextdigit;
}

string generate_next_number(string s, static string mapdigit = "2235577222") {
  if (s.length() == 0)
    return "2";
  char digit = s.back();
  char nextdigit = mapdigit[digit - '0'];
  s.pop_back();
  if (nextdigit < digit)
    s = generate_next_number(s);
  s.push_back(nextdigit);
  return s;
}

// https://stackoverflow.com/questions/74134479/given-a-number-find-the-next-higher-number-containing-only-certain-digits/74388175#74388175
