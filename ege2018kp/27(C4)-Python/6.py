# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-6. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На вход программы подаются сведения о результатах соревнований по
школьному многоборью. Многоборье состоит из соревнований по четырем
видам спорта, участие в каждом из которых оценивается баллами от 0 до
10 (0 баллов получает ученик, не принимавший участия в соревнованиях
по данному виду спорта). Победители определяются по наибольшей сумме 
набранных баллов. Известно, что общее количество участников
соревнований не превосходит 100.
В первой строке вводится количество учеников, принимавших участие в
соревнованиях, N. Далее следуют N строк, имеющих следующий формат: 
    <Фамилия> <Имя> <Баллы>
Здесь <Фамилия> – строка, состоящая не более чем из 20 символов; <Имя>
– строка, состоящая не более чем из 15 символов; <Баллы> - строка, 
содержащая четыре целых числа, разделенных пробелом, соответствующих
баллам, полученным на соревнованиях по каждому из четырех видов
спорта. При этом <Фамилия> и <Имя>, <Имя> и <Баллы> разделены одним
пробелом. Примеры входных строк:    
    Иванова Мария 5 8 6 3
    Петров Сергей 9 9 5 7
Напишите программу, которая будет выводить на экран фамилии и имена
трех лучших участников многоборья. Если среди остальных участников 
есть ученики, набравшие то же количество баллов, что и один из трех
лучших, то их фамилии и имена также следует вывести. При этом имена и
фамилии можно выводить в произвольном порядке.
"""
import sys, codecs
save_stdin = sys.stdin
sys.stdin = codecs.open("in/6.in", "r", "utf-8")

N = int(input())
res = {}
for i in range(N):
    fam, name, b1, b2, b3, b4 = input().split()
    ball = [int(x) for x in [b1, b2, b3, b4]]
    res[fam+" "+name] = sum(ball)
res = sorted(res.items(), key = lambda x: x[1], reverse = True)
best = [x[0] for x in res if x[1] >= res[2][1]]
for x in best: 
    print(x)

sys.stdin = save_stdin