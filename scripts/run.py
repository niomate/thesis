import glob
import queue
import sys
import time
from multiprocessing import Pool, Process, Queue, cpu_count, current_process
from os.path import abspath, dirname, join
from subprocess import DEVNULL, run


def log_thread(*args):
    print(f'[{current_process().name}]', *args)


def inpaint(inim, mask, *args, outim=None):
    log_thread(f'Inpainting {inim}')
    if outim is None:
        outim = inpaint_name(inim)
    cmd = ['build/inpainting', inim, mask, *args, '-o', outim]
    print(' '.join(cmd))
    run(cmd)


def corners(inim, *args, outim=None):
    ''' Pass all arguments that build/corners accepts '''
    log_thread(f'Corner detection {inim}')
    if outim is None:
        outim = mask_name(inim)
    cmd = ['build/corners', inim, '-o', outim, *args]
    print(' '.join(cmd))
    run(cmd)  # stdout=DEVNULL)
    return outim


def name_from_path(path):
    return path.rsplit('/', 1)[1].split('.')[0]


def mask_name(path):
    n = name_from_path(path)
    return f'output/masks/{n}_mask.pgm'


def inpaint_name(path):
    n = name_from_path(path)
    return f'output/inpainted/{n}_inpainted.pgm'


def worker(todo, mode):
    while True:
        if todo.empty():
            log_thread('Done')
            break
        try:
            image = todo.get(False)
            if mode == 'full':
                mask = corners(image, *sys.argv[2:])
                inpaint(image, mask)
            elif mode == 'corners':
                corners(image, *sys.argv[2:])
            elif mode == 'mask':
                corners(image, *sys.argv[2:])
        except queue.Empty:
            continue


if __name__ == '__main__':

    if len(sys.argv) == 1:
        mode = 'full'
    else:
        mode = sys.argv[1]

    inimages = glob.iglob('images/binary/**/*.pgm', recursive=True)
    todo = Queue()

    for im in inimages:
        todo.put(im)

    jobs = []
    for i in range(6):
        p = Process(target=worker, args=(todo, mode))
        p.name = f'Worker-{i}'
        p.start()

    for p in jobs:
        p.join()
