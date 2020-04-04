sequence = list(map(int, input().split()))
max_num, max_len, cur_len = 0, 0, (1 if sequence[0] >= sequence[1] else 0)
for i in range(len(sequence) - 1):
    if cur_len > 0 and sequence[i] == sequence[i - 1] and sequence[i] >= sequence[i + 1]:
        cur_len += 1
    elif sequence[i - 1] <= sequence[i] >= sequence[i + 1]:
        cur_len = 1
    else:
        cur_len = 0
    max_len, max_num = (cur_len, sequence[i]) if cur_len > max_len else (max_len, max_num)

if cur_len > 0 and sequence[-1] == sequence[-2]:
    cur_len += 1
elif sequence[i - 1] <= sequence[i] >= sequence[i + 1]:
    cur_len = 1
max_len, max_num = (cur_len, sequence[-1]) if cur_len > max_len else (max_len, max_num)

print(max_num)

max_num = cur_len = 0
with open("t.txt", "r", encoding="utf-8") as f:
    prev1 = int(f.readline())
    prev2 = int(f.readline())
    cur = f.readline()
    cur_len = 1 if prev1 >= prev2 else 0
    while cur:
        if cur_len > 0 and prev2 == prev1 and prev2 >= int(cur):
            cur_len += 1
        elif prev1 <= prev2 <= int(cur):
            cur_len = 1
        else:
            cur_len = 0
        max_len, max_num = (cur_len, prev2) if cur_len > max_len else (max_len, max_num)
        prev1, prev2, cur = prev2, int(cur), int(f.readline())

if cur_len > 0 and prev2 == prev1:
    cur_len += 1
elif sequence[i - 1] <= sequence[i] >= sequence[i + 1]:
    cur_len = 1
max_len, max_num = (cur_len, sequence[-1]) if cur_len > max_len else (max_len, max_num)

print(max_num)
