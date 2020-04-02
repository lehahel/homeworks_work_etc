sequence = list(map(int, input().split()))
x = int(input())
print(sum(list(sequence[i] * (x ** i) for i in range(len(sequence)))))
print(sum(list(sequence[i] * i * (x ** (i - 1)) for i in range(1, len(sequence)))))

