# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-54. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2014
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На вход программы подаются результаты измерений активности 
космических частиц, выполняемых прибором с интервалом 1 минуту. 
Все данные – целые числа. Требуется найти наибольшую сумму двух 
результатов измерений, выполненных с интервалом не менее, чем 
в 7 минут. 
Описание входных данных
  В первой строке вводится одно целое положительное число – 
  количество измерений N, которое может быть очень велико. 
  Гарантируется, что N > 7. Каждая из следующих N строк содержит 
  по одному целому числу – результат очередного измерения.
Описание выходных данных
  Программа должна вывести одно число наибольшую сумму двух 
  результатов измерений, выполненных с интервалом не менее, 
  чем в 7 минут.
Пример входных данных:
  10
  1
  2
  3
  4
  5
  6
  7
  8
  9
  10
Пример выходных данных для приведённого выше примера входных данных:
  13
"""
import sys
save_stdin = sys.stdin
sys.stdin = open("in/54.in")

N = int(input())
K = 7
Buf = [0]*K
for i in range(N):
  elem = int(input())
  next = Buf[i % K] 
  Buf[i % K] = elem
  if i == K:
    maxPrev = next
    maxSum = next + elem
  if i > K:
    maxPrev = max(maxPrev, next)
    maxSum = max(maxSum, maxPrev+elem)

print(maxSum)

sys.stdin = save_stdin