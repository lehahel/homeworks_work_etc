# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-7. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
В некотором вузе абитуриенты проходят предварительное тестирование, по
результатам которого могут быть допущены к сдаче вступительных
экзаменов в первом потоке. Тестирование проводится по двум предметам,
по каждому предмету абитуриент может набрать от 0 до 100 баллов. При
этом к сдаче экзаменов в первом потоке допускаются абитуриенты,
набравшие по результатам тестирования не менее 30 баллов по каждому из
двух предметов. На вход программы подаются сведения о результатах
предварительного тестирования. Известно, что общее количество
участников тестирования не превосходит 500. 
В первой строке вводится количество абитуриентов, принимавших участие
в тестировании, N. Далее следуют N строк, имеющих следующий формат: 
    <Фамилия> <Имя> <Баллы> 
Здесь <Фамилия> – строка, состоящая не более чем из 20 символов; <Имя>
– строка, состоящая не более чем из 15 символов; <Баллы> – строка,
содержащая два целых числа, разделенных пробелом, соответствующих
баллам, полученным на тестировании по каждому из двух предметов. При
этом <Фамилия> и <Имя>, <Имя> и <Баллы> разделены одним пробелом.
Примеры входных строк:  
    Ветров Роман 68 59
    Анисимова Екатерина 64 88 
Напишите программу, которая будет выводить на экран фамилии и имена
абитуриентов, потерпевших неудачу, то есть не допущенных к сдаче
экзаменов в первом потоке. При этом фамилии должны выводиться в
алфавитном порядке.
"""
import sys, codecs
save_stdin = sys.stdin
sys.stdin = codecs.open("in/7.in", "r", "utf-8")

N = int(input())
failed = []
for i in range(N):
    fam, name, b1, b2 = input().split()
    ball = [int(x) for x in [b1, b2]]
    if ball[0] < 30 or ball[1] < 30:
        failed.append(fam+" "+name)
for x in sorted(failed): 
    print(x)

sys.stdin = save_stdin