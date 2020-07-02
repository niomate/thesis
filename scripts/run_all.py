import glob
import queue
import sys
from multiprocessing import Process, Queue, current_process
from subprocess import run


def log_thread(*args):
    print(f'[{current_process().name}]', *args)


def str_list(l):
    return list(map(str, l))


def inpaint(inim, mask, *args, outim=None):
    log_thread(f'Inpainting {inim}')
    args = str_list(args)
    if outim is None:
        outim = inpaint_name(inim)
    cmd = ['./inpaint', inim, mask, *args, '-o', outim]
    print(' '.join(cmd))
    run(cmd)


def mask(inim, *args, outim=None):
    ''' Pass all arguments that corner accepts '''
    log_thread(f'Mask computation {inim}')
    args = str_list(args)
    if outim is None:
        outim = mask_name(inim)
    cmd = ['./corner', inim, '-M', '-o', outim, *args]
    print(' '.join(cmd))
    run(cmd)
    return outim


def corners(inim, *args, outim=None):
    ''' Pass all arguments that corner accepts '''
    log_thread(f'Corner detection {inim}')
    args = str_list(args)
    if outim is None:
        outim = mask_name(inim)
    cmd = ['./corner', inim, '-o', outim, *args]
    print(' '.join(cmd))
    run(cmd)
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
                mask(image, *sys.argv[2:])
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
