﻿# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-34. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На вход программе подается текстовый файл, в первой строке которого
записано квадратное уравнение, причем используются обозначения:
    x^2 обозначается как «a»
    x обозначается как «b»
Например, уравнение 2*x^2-4*x-6=0  запишется в виде строки
    2a-4b-6
Гарантируется, что уравнение имеет «хороший» вид, все его коэффициенты
определены и корни вещественные. Напишите эффективную по времени работы 
и по используемой памяти программу (укажите используемую версию языка
программирования, например, Borland Pascal 7.0), которая выводит на экран
корни уравнения. Для приведенного входного файла программа
должна вывести
    -1
    3 
"""
import sys
save_stdin = sys.stdin
sys.stdin = open("in/34.in")

s = input() + '+'

import math
coef = [0, 0, 0]
for i in range(3):
    k = 1
    while not s[k] in "+-": k += 1
    if s[k-1] == 'a': 
        coef[2] = int(s[0:k-1])
    elif s[k-1] == 'b': 
        coef[1] = int(s[0:k-1])
    else: coef[0] = int(s[0:k])
    s = s[k:]
    
a, b, c = coef[2], coef[1], coef[0]
D = b*b - 4*a*c
print("%.3f" % ( (-b + math.sqrt(D)) / (2*a)) )
print("%.3f" % ( (-b - math.sqrt(D)) / (2*a)) )

sys.stdin = save_stdin