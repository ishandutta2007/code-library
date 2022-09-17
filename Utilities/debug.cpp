
template <class T> void debug(std::initializer_list<T> list) {
  #ifndef ONLINE_JUDGE
  for (auto elem : list)
    cout << elem << " ";
  cout << endl;
  #endif
}