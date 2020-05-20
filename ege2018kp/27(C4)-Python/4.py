# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-4. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На вход программы подаются фамилии и имена учеников. Известно, что
общее количество учеников не превосходит 100. В первой строке вводится
количество учеников, принимавших участие в соревнованиях, N. Далее
следуют N строк, имеющих следующий формат: 
    <Фамилия> <Имя> 
Здесь <Фамилия> – строка, состоящая не более чем из 20 символов; <Имя>
– строка, состоящая не более чем из 15 символов. При этом <Фамилия> и
<Имя> разделены одним пробелом. Примеры входных строк:  
    Иванова Мария
    Петров Сергей 
Требуется написать программу, которая формирует и печатает уникальный
логин для каждого ученика по следующему правилу: если фамилия
встречается первый раз, то логин – это данная фамилия, если фамилия
встречается второй раз, то логин – это фамилия, в конец которой
приписывается число 2 и т.д. Например, для входной последовательности 
    Иванова Мария 
    Петров Сергей 
    Бойцова Екатерина 
    Петров Иван 
    Иванова Наташа
будут сформированы следующие логины: 
    Иванова 
    Петров
    Бойцова
    Петров2 
    Иванова2
"""
import sys, codecs
save_stdin = sys.stdin
sys.stdin = codecs.open("in/4.in", "r", "utf-8")

N = int(input())
login = {}
for i in range(N):
    fam, name = input().split()
    if login.get(fam) == None:
        login[fam] = 1
        print(fam)
    else:
        login[fam] += 1
        print("%s%d" % (fam, login[fam]))

sys.stdin = save_stdin