import cipher
import file_work
import trainer
import argparser


args = argparser.parse_arguments()
if args.func == 'encode':
    text = file_work.read_text(args.input_file)
    result = cipher.encode(args.cipher, int(args.key), text)
    file_work.write_text(result, args.output_file)

elif args.func == 'decode':
    text = file_work.read_text(args.input_file)
    result = cipher.decode(args.cipher, int(args.key), text)
    file_work.write_text(result, args.output_file)

elif args.func == 'train':
    text = file_work.read_text(args.text_file)
    trainer.train(text, args.model_file)

else:
    text = file_work.read_text(args.input_file)
    result = trainer.caesar_hack(text, args.model_file)
    file_work.write_text(result, args.output_file)
