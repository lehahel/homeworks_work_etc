sequence = list(map(int, input().split()))
print(sum(list(sequence.count(x) for x in range(1, 6))))

