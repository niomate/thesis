import os
import sys
import progressbar
from contextlib import contextmanager
from ctypes import cdll, c_char_p
from os.path import join, abspath, dirname


@contextmanager
def stdout_redirected(to=os.devnull):
    fd = sys.stdout.fileno()

    def __redirect_stdout(to):
        sys.stdout.close()
        os.dup2(to.fileno(), fd)
        sys.stdout = os.fdopen(fd, 'w')

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'w') as file:
            __redirect_stdout(to=file)
        try:
            yield
        finally:
            __redirect_stdout(to=old_stdout)


here = dirname(__file__)
images = abspath(join(here, '../images/binary/angles/'))
output = abspath(join(here, '../output/{}-out.pgm')).format
lib = cdll.LoadLibrary(abspath(join(here, '../lib/libamss.so')))
amss = lib.test_detection
amss.argtypes = [c_char_p, c_char_p]

nimg = len(os.listdir(images))
bar = progressbar.ProgressBar(max_value=nimg)

i = 0
for image in os.listdir(images):
    name = image.rsplit('.', 1)[0]
    out = output(name)
    inpath = join(images, image)
    # print(f'Testing {name}...')
    with stdout_redirected():
        amss(inpath.encode('utf-8'), out.encode('utf-8'))
    bar.update(i)
    i += 1
