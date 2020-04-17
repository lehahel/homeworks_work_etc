import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='func', help='encode/decode/train/hack')

    encode_parser = subparsers.add_parser("encode")
    decode_parser = subparsers.add_parser("decode")
    train_parser = subparsers.add_parser("train")
    hack_parser = subparsers.add_parser("hack")

    encode_parser.add_argument("--cipher", choices=['caesar', 'vigenere'], dest='cipher', help="type of cipher you use")
    decode_parser.add_argument("--cipher", choices=['caesar', 'vigenere'], dest='cipher', help="type of cipher you use")
    encode_parser.add_argument("--key", dest='key', help="key of cipher: number for Caesar, string for Vigenere")
    decode_parser.add_argument("--key", dest='key', help="key of cipher: number for Caesar, string for Vigenere")

    encode_parser.add_argument("--input-file", dest='input_file', default=None, help="text file to encode/decode")
    decode_parser.add_argument("--input-file", dest='input_file', default=None, help="text file to encode/decode")
    hack_parser.add_argument("--input-file", dest='input_file', default=None, help="text file to encode/decode")

    encode_parser.add_argument("--output-file", dest='output_file', default=None,
                               help="output file to encode/decode into")
    decode_parser.add_argument("--output-file", dest='output_file', default=None,
                               help="output file to encode/decode into")
    hack_parser.add_argument("--output-file", dest='output_file', default=None,
                             help="output file to encode/decode into")

    train_parser.add_argument("--text-file", dest='text_file', default=None, help="text file to train")

    train_parser.add_argument("--model-file", dest='model_file', help="file with statistic to hack Caesar cipher")
    hack_parser.add_argument("--model-file", dest='model_file', help="file with statistic to hack Caesar cipher")

    args = parser.parse_args()
    return args
