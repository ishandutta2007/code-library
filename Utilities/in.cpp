//equivalent of python in
bool in(int value, vector<int>container) {
  auto it = std::find(container.begin(), container.end(), value);
  return (it != container.end());
}
