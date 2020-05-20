# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-67. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2017
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
67) На вход программы поступает последовательность 
из N натуральных чисел. Требуется определить, какая 
цифра чаще всего встречается в десятичной записи 
этих чисел. Если таких цифр несколько, необходимо 
вывести их все в порядке убывания – от большей к 
меньшей.
Входные данные:
  На вход программе подаётся натуральное число 
  N (N <= 1000), а затем N натуральных чисел, каждое 
  из которых не превышает 10000. 
Пример входных данных: 
3 
13
214
32
Выходные данные:
  Программа должна вывести цифры, которые встречаются 
  в последовательности наибольшее число раз, в 
  порядке убывания. 
Пример выходных данных для приведённого примера входных данных: 
3 2 1
(цифры 3, 2 и 1 встречаются по 2 раза).
 
"""
import sys
save_stdin = sys.stdin
sys.stdin = open("in/67.in")

count = [0]*10
def Digits(x):
  global count
  while x > 0:
    d = x % 10
    count[d] += 1
    x //= 10

N = int(input())
for i in range(N):
  x = int(input())
  Digits(x)
  
M = max(count)
arrD = [x for x in range(9,-1,-1) 
          if count[x] == M]
for d in arrD:
  print(d, end=" ")

sys.stdin = save_stdin