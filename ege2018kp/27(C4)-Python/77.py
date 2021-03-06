# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-77. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2017
-------------------------------------------------
77) (Д.В. Богданов) Дан набор из N натуральных чисел. 
Необходимо определить количество пар элементов (ai, aj) 
этого набора, в которых 1 ≤ i < j ≤ N и сумма 
элементов кратна 12. Напишите эффективную по времени и 
по памяти программу для решения этой задачи. 
Описание входных и выходных данных 
В первой строке входных данных задаётся количество чисел 
N (1 ≤ N ≤ 10000). В каждой из последующих N строк записано 
одно натуральное число, не превышающее 1000.
Пример входных данных:
  5
  7
  5
  6
 12
 24
Пример выходных данных для приведённого выше примера 
входных данных:
  2
В приведённом наборе из 5 чисел имеются две пары 
(7, 5), (12, 24), сумма элементов которых кратна 12.
"""
import sys
save_stdin = sys.stdin
sys.stdin = open("in/77.in")

rem = [0]*12
N = int(input())
for i in range(N):
  x = int(input())
  rem[x%12] += 1

count = (rem[0]*(rem[0]-1) +
         rem[6]*(rem[6]-1)) // 2
for i in range(1, 6):
  count += rem[i]*rem[12-i];

print(count)

sys.stdin = save_stdin