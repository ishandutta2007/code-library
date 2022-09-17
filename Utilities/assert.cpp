
#define ASSERT(condition)                                                      \
  {                                                                            \
    if (!(condition)) {                                                        \
      cout << "ASSERT FAILED: " << #condition << " @ " << __FILE__ << " ("     \
           << __LINE__ << ")" << std::endl;                                    \
    }                                                                          \
    exit(0);                                                                   \
  }
