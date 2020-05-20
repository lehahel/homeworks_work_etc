# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-8. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На вход программе подаются сведения о телефонах всех сотрудников
некоторого учреждения. В первой строке сообщается количество
сотрудников N, каждая из следующих N строк имеет следующий формат: 
    <Фамилия> <Инициалы> <телефон>
где <Фамилия> – строка, состоящая не более чем из 20 символов,
<Инициалы> - строка, состоящая не более чем из 4-х символов (буква,
точка, буква, точка), <телефон> – семизначный номер, 3-я и 4, я, а
также 5-я и 6-я цифры которого разделены символом «–». <Фамилия> и
<Инициалы>, а также <Инициалы> и <телефон> разделены одним пробелом.
Пример входной строки:
    Иванов П.С. 555-66-77
Сотрудники одного подразделения имеют один и тот же номер телефона.
Номера телефонов в учреждении отличаются только двумя последними
цифрами. Требуется написать как можно более эффективную программу,
которая будет выводить на экран информацию, сколько в среднем
сотрудников работает в одном подразделении данного учреждения.
"""
import sys, codecs
save_stdin = sys.stdin
sys.stdin = codecs.open("in/8.in", "r", "utf-8")

N = int(input())
podr = {}
for i in range(N):
    fam, io, tel = input().split()
    x, x, code = tel.split('-')
    podr[code] = podr.get(code,0) + 1
numPeople = [x[1] for x in podr.items()]
print("%.3f" % (sum(numPeople)/len(numPeople)) )

sys.stdin = save_stdin