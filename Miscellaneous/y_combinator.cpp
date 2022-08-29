template <typename Func> class YCombinatorResult {
public:
  template <typename T>
  explicit YCombinatorResult(T &&func)
      : func_(std::forward<T>(func)) {}

  template <typename... Args> decltype(auto) operator()(Args &&... args) {
    return func_(std::ref(*this), std::forward<Args>(args)...);
  }

private:
  Func func_;
};

template <typename Func> decltype(auto) y_combinator(Func &&func) {
  return YCombinatorResult<std::decay_t<Func>>(std::forward<Func>(func));
}

void example() {
  auto fact = y_combinator(
      [&](auto self, int x) -> int { return x == 0 ? 1 : self(x - 1) * x; });
  std::cerr << fact(10) << std::endl;
}

// https://github.com/ssstare/icpc/tree/fec306e49c04d516d1a2ef1ae0c8b20e764c543f/snippets