import numpy as np
import subprocess
import os
import sys
from PIL import Image

here = os.path.dirname(__file__)
harris = os.path.join(here, '../harris-corner-detector_1.0/bin/harris_corner_detector')
results = os.path.join(here, '../.tmp')


def get_corners(img):
    with open(os.path.join(results, f'{img}.corners')) as f:
        lines = f.readlines()
    corners = [[float(x) for x in l.split()[:-1]] for l in lines[1:]]
    return corners


def harris_detector(img, verbose=True):
    name = img.rsplit('/', maxsplit=1)[1].split('.')[0]
    cmd = [harris, img, '-f', os.path.join(here, f'../.tmp/{name}.corners')]
    print(f'Running {" ".join(cmd)}')
    if verbose:
        cmd.append('-v')
    subprocess.run(cmd)
    corners = get_corners(name)
    m = mask(img, corners)
    m.save(os.path.join(here, f'../.tmp/{name}_mask.pgm'))


def clamp(i, hi, lo):
    if i > hi:
        return hi
    if i < lo:
        return lo
    return i


def mask(image_path, corners, radius=4):
    im = Image.open(image_path)
    image = np.array(im).T
    w, h = image.shape
    new = np.zeros((w, h))
    for y, x in corners:
        x = int(np.round(x))
        y = int(np.round(y))
        for i in range(x-radius-1, x+radius+1):
            for j in range(y-radius-1, y+radius+1):
                if ((x-i)**2+(y-j)**2) <= radius**2:
                    nx = clamp(i, w, 0)
                    ny = clamp(j, h, 0)
                    new[nx, ny] = 255
    m = Image.fromarray(new.astype(np.uint8), mode='L')
    m.show()
    return m


if __name__ == '__main__':
    harris_detector(sys.argv[1])
