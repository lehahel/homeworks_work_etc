import file_work
import cipher
import json


alphabet_size = 26


def hack(text, cipher_type, model_file):
    if cipher_type == 'caesar':
        return caesar_hack(text, model_file)
    else:
        raise NameError('Wrong cipher name')


def get_stat(text):
    stat = [0 for i in range(alphabet_size)]
    full_sum = 0
    for x in text:
        if ord('a') <= ord(x) <= ord('z'):
            stat[ord(x) - ord('a')] += 1
            full_sum += 1
    stat = [x / full_sum for x in stat]
    return stat


def find_model_distance(stat1, stat2):
    if len(stat1) != alphabet_size or len(stat2) != alphabet_size:
        raise NameError('Wrong stat')
    distance = 0
    for i in range(alphabet_size):
        distance += abs(stat1[i] - stat2[i])
    return distance


def shift(stat, value):
    if len(stat) != alphabet_size:
        raise NameError('Wrong stat')
    new_stat = [0 for i in range(alphabet_size)]
    first_el = stat[0]
    for i in range(alphabet_size):
        new_stat[i] = stat[(i + value) % alphabet_size] if (i + value) % alphabet_size != 0 else first_el
    return new_stat


def caesar_hack(text, model_file):
    my_stat = get_stat(text)
    with open(model_file, "r") as f:
        true_stat = json.load(f)
    best_k = 0
    best_distance = -1
    for k in range(1, alphabet_size):
        shifted_stat = shift(my_stat, k)
        distance = find_model_distance(shifted_stat, true_stat)
        if distance < best_distance or best_distance == -1:
            best_distance = distance
            best_k = k
    return cipher.decode('caesar', best_k, text)


def read_stat_file(filename):
    if not isinstance(filename, str):
        raise TypeError('filename should be str')
    if not ('.txt' in filename):
        raise NameError('file should be .txt')
    stat = []
    file = open(filename)
    for i in range(alphabet_size):
        stat.append(float(file.readline().replace('\n', '')))
    file.close()
    return stat


def train(text, model_file):
    stat = get_stat(text)
    with open(model_file, "w") as f:
        json.dump(stat, f)
