﻿# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-42. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На ускорителе для большого числа частиц производятся замеры скорости
каждой из них. Чтобы в документации качественно отличать одну серию от
другой, каждую серию решили характеризовать числом, равным минимальному 
произведению из всех произведений пар скоростей различных частиц. Вам
предлагается написать эффективную, в том числе по используемой памяти,
программу (укажите используемую версию языка программирования, например
Borland Pascal 7.0), которая будет обрабатывать результаты эксперимента,
находя искомую величину. В нашей модели скорость частицы - это величина,
которая может принимать как положительные, так и отрицательные значения.
Следует учитывать, что частиц, скорость которых измерена, может быть
очень много, но не может быть меньше двух. 
Перед текстом задачи кратко опишите используемый вами алгоритм решения 
задачи. 
На вход программе в первой строке подается количество частиц N. В каждой
из последующих N строк записано одно целое число со знаком (плюс или 
минус), по абсолютной величине не превосходящее 10000.
Пример входных данных: 
    5 
    +123 
    +2000 
    +10 
    +3716 
    +10 
Программа должна вывести одно число - минимальное произведение из всех
произведений пар скоростей различных частиц.
Пример выходных данных для приведенного выше примера входных данных: 
    100
"""
import sys
save_stdin = sys.stdin
sys.stdin = open("in/42.in")

N = int(input())
v1 = int(input())
v2 = int(input())
maX, max2 = max(v1,v2), min(v1,v2)
miN, min2 = min(v1,v2), max(v1,v2)
for i in range(N-2):
    v = int(input())
    if v < miN:
        min2 = miN
        miN = v
    elif v < min2: min2 = v        
    if v > maX:
        max2 = maX 
        maX = v
    elif v > max2: max2 = v

print(min([miN*min2, maX*miN, maX*max2]))

sys.stdin = save_stdin