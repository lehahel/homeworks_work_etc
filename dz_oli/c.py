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
