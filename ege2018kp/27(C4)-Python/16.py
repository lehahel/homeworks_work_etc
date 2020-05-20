# -*- coding: utf-8 -*-
"""
-------------------------------------------------
Задача C4-16. 
Решение на языке Python 3.
 Автор: Константин Поляков, 2013
E-mail: kpolyakov@mail.ru
   Web: kpolyakov.spb.ru
-------------------------------------------------
На вход программе подаются сведения о пассажирах, сдавших свой багаж в
камеру хранения. В первой строке задано текущее время: через двоеточие
два целых числа, соответствующие часам (от 00 до 21, ровно 2 символа)
и минутам (от 00 до 59, ровно 2 символа). Во второй строке задается
количество пассажиров N, которое не меньше 10, но не превосходит 1000.
В каждой из последующих N строк находится информация о пассажирах в
следующем формате: 
  <Фамилия> <Время освобождения ячейки> 
где <Фамилия> – строка, состоящая не более, чем из 20 символов без
пробелов, <Время освобождения ячейки> – через двоеточие два целых
числа, соответствующие часам (от 00 до 21, ровно 2 символа) и минутам
(от 00 до 59, ровно 2 символа). <Фамилия> и <Время освобождения 
ячейки> разделены ровно одним пробелом. Пример входных строк: 
  10:00 
  3
  Иванов 12:00
  Петров 10:12
  Сидоров 12:12 
Программа должна выводить список пассажиров, которые в ближайшие 2 
часа должны освободить ячейки. Список должен быть отсортирован в
хронологическом порядке освобождения ячеек. В данном случае программа
должна вывести 
  Петров
  Иванов
"""
import sys, codecs
save_stdin = sys.stdin
sys.stdin = codecs.open("in/16.in", "r", "utf-8")

def time2Int ( t ):
    h, m = t.split(':')
    return 60*int(h) + int(m)

curTime = time2Int(input())
N = int(input())
free = []
for i in range(N):
    fam, t = input().split()
    tInt = time2Int(t)
    if tInt <= curTime + 120:
        free.append( (fam, t) )
for x in sorted(free, key = lambda x: x[1]):
    print(x[0])

sys.stdin = save_stdin