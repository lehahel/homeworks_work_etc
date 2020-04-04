result = 0
deriv = 0
with open("t.txt", "r", encoding="utf-8") as f:
    degree = 1
    x = int(f.readline())
    result += int(f.readline())
    cur = f.readline()
    while cur:
        result += int(cur) * (x ** degree)
        deriv += int(cur) * degree * (x ** (degree - 1))
        degree += 1
        cur = f.readline()

with open("out.txt", "w", encoding="utf-8") as f:
    f.write("{} {}".format(str(result), str(deriv)))
