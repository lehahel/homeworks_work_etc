# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-65. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2016
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
64) Задание Б. Имеется набор данных, состоящий из пар 
положительных целых чисел. Необходимо выбрать из каждой 
пары ровно одно число так, чтобы сумма всех выбранных чисел 
не делилась на 4 и при этом была максимально возможной. Если 
получить требуемую сумму невозможно, в качестве ответа нужно 
выдать 0.
Программа считается эффективной по времени, если время работы 
программы пропорционально количеству пар чисел N, т.е. при 
увеличении N в k раз время работы программы должно увеличиваться 
не более чем в k раз. Программа считается эффективной по памяти, 
если размер памяти, использованной в программе для хранения данных, 
не зависит от числа N и не превышает 1 килобайта.
На вход программе в первой строке подаётся количество 
пар N (1 <= N <= 100000). Каждая из следующих N строк 
содержит два натуральных числа, не превышающих 10 000. 
Пример входных данных: 
6
1 3 
5 12 
6 8 
5 4 
3 3 
1 1
Пример выходных данных для приведённых выше примеров входных данных: 
31

"""
import sys
save_stdin = sys.stdin
sys.stdin = open("in/65.in")

N = int(input())

sum = 0  
delta = 10001
for i in range(N):
  x1, x2 = map(int, input().split())
  sum += max(x1,x2)
  d = abs(x1-x2)
  if d % 4 != 0 and d < delta:
     delta = d
if sum % 4 == 0:
  if delta == 10001:
        sum = 0  
  else: sum = sum - delta
print(sum);

sys.stdin = save_stdin