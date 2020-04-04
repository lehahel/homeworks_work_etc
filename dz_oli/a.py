res = 0
with open("t.txt", "r", encoding="utf-8") as f:
    x = f.readline()
    while x:
        if 1 <= int(x) <= 5:
            res += 1
        x = f.readline()

with open("out.txt", "w", encoding="utf-8") as f:
    f.write(str(res))
