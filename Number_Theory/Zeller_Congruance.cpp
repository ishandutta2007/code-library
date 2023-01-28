#include <iostream>
#include <cmath>
using namespace std;

struct date {
  int day;
  int month;
  long long int year;
  int day_of_week;
  date(int a, int b, long long int c) {
    day = a;
    month = b;
    year = c;
    day_of_week = (*this).zeller_congruence();
  }
  void increment_by_one() {
    day_of_week = (day_of_week + 1) % 7;
    bool type31 = (month == 1) || (month == 3) || (month == 5) ||
                  (month == 7) || (month == 8) || (month == 10) ||
                  (month == 12);
    bool type30 = (month == 4) || (month == 6) || (month == 6) ||
                  (month == 9) || (month == 11);
    bool typefeb = month == 2;
    bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    bool month_end = (type31 && day == 31) || (type30 && day == 30) ||
                     (typefeb && day == 28 && !leap) ||
                     (typefeb && day == 29 && leap);
    bool year_end = (month_end && month == 12);

    if (year_end) {
      day = 1;
      month = 1;
      year++;
    } else if (month_end) {
      day = 1;
      month++;
    } else
      day++;
  }

  int zeller_congruence() {
    int q = day;
    int m = month;
    long long int y = year;
    if (m == 1) {
      m = 13;
      y--;
    }
    if (m == 2) {
      m = 14;
      y--;
    }

    int k = y % 100;
    long long int j = y / 100;
    return ((q + (13 * (m + 1)) / 5 + k + k / 4 + j / 4 + 5 * j) % 7);
  }

  void print_date() {
    cout << day << " " << month << " " << year << " ->" << day_of_week << endl;
  }
};

bool isDateGreaterthan(date a, date b) {
  if (a.year > b.year)
    return true;
  else if (a.year < b.year)
    return false;
  else {
    if (a.month > b.month)
      return true;
    else if (a.month < b.month)
      return false;
    else {
      if (a.day > b.day)
        return true;
      else if (a.day < b.day)
        return false;
      else
        return false;
    }
  }
}

int countSundays(date dt1, date dt2) {
  int count = 0;
  while ((dt1.day != dt2.day) || (dt1.month != dt2.month) ||
         (dt1.year != dt2.year)) {
    if (dt1.day == 1 && dt1.day_of_week == 1) {
      count++;
      // dt1.print_date();
    }
    dt1.increment_by_one();
  }
  if (dt1.day == 1 && dt1.day_of_week == 1)
    count++;
  return count;
}

int main() {
  int t;
  cin >> t;
  for (int a0 = 0; a0 < t; a0++) {
    int d1, m1, d2, m2;
    long long int y1, y2;
    cin >> y1 >> m1 >> d1;
    cin >> y2 >> m2 >> d2;
    date dt1 = date(d1, m1, y1);
    date dt2 = date(d2, m2, y2);
    // cout << dt1.zeller_congruence()<<endl;
    // cout << dt2.zeller_congruence()<<endl;
    if (isDateGreaterthan(dt1, dt2)) {
      cout << countSundays(dt2, dt1) << endl;
    } else
      cout << countSundays(dt1, dt2) << endl;
  }
}