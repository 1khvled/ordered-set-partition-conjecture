p = 'hunt.py'
b = open(p, 'rb').read()
idx = b.rfind(b'if __name__')
b = b[:idx] + b'if __name__ == "__main__":\n    main()\n'
open(p, 'wb').write(b)
print('done')
