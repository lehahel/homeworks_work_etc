alphabet_size = 26


def encode(cipher_type, key, text):
    if cipher_type == 'caesar':
        return caesar_shift(key, text, 'encode')
    elif cipher_type == 'vigenere':
        return vigenere_shift(key, text, 'decode')


def decode(cipher_type, key, text):
    if cipher_type == 'caesar':
        return caesar_shift(key, text, 'decode')
    elif cipher_type == 'vigenere':
        return vigenere_shift(key, text, 'decode')


def caesar_shift(key, text, mode):
    shift = key if mode == 'encode' else -key
    new_text = ""
    for ch in text:
        if ord('A') <= ord(ch) <= ord('Z'):
            new_text += chr((ord(ch) - ord('A') + alphabet_size + shift) % alphabet_size + ord('A'))
        elif ord('a') <= ord(ch) <= ord('z'):
            new_text += chr((ord(ch) - ord('a') + alphabet_size + shift) % alphabet_size + ord('a'))
        else:
            new_text += ch
    return new_text


def vigenere_shift(key, text, mode):
    new_text = ""
    additions = list(ord(x) - ord('a') for x in key) if mode == 'encode' else list(-ord(x) + ord('a') for x in key)
    for i in range(len(text)):
        if ord('A') <= ord(text[i]) <= ord('Z'):
            new_text += chr((ord(text[i]) - ord('A') + alphabet_size + additions[i % len(additions)])
                            % alphabet_size + ord('A'))
        elif ord('a') <= ord(text[i]) <= ord('z'):
            new_text += chr((ord(text[i]) - ord('a') + alphabet_size + additions[i % len(additions)])
                            % alphabet_size + ord('a'))
        else:
            text += text[i]
