import string

digs = string.digits + string.ascii_letters


def int2base(x, base):
    if x < 0:
        sign = -1
    elif x == 0:
        return digs[0]
    else:
        sign = 1

    x *= sign
    digits = []

    while x:
        digits.append(digs[x % base])
        x = x // base

    if sign < 0:
        digits.append("-")
    digits.reverse()
    return "".join(digits)


# # Base 10
# def isSelfDescribing10(n):
#     s = str(n)
#     return all(s.count(str(i)) == int(ch) for i, ch in enumerate(s))


# for i in range(1, 1000001):
#     if isSelfDescribing10(i):
#         print(i, isSelfDescribing10(i))

# ans = [
#     (x, isSelfDescribing10(x))
#     for x in (1210, 2020, 21200, 3211000, 42101000, 521001000, 6210001000)
# ]
# print(ans)


# # Base 4
# def isSelfDescribing4(n):
#     s = str(int2base(n, 4))
#     return all(s.count(str(i)) == int(ch) for i, ch in enumerate(s))


# for i in range(1, 1000001):
#     if isSelfDescribing4(i):
#         print(i, isSelfDescribing4(i))


# Any base
def isSelfDescribing_any_base(n):
    for b in range(2, 11):
        s = str(int2base(n, b))
        if len(s) != b:
            continue
        is_sd = all(s.count(str(i)) == int(ch) for i, ch in enumerate(s))
        if is_sd:
            return is_sd, b
    return False, False


for i in range(1, 1000001):
    is_sd_any, base = isSelfDescribing_any_base(i)
    if is_sd_any:
        print(i, is_sd_any, "in base=(", base, ")")
