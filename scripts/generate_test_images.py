from PIL import Image
from PIL import ImageDraw
import os
import sys
import numpy as np

if __name__ == '__main__':
    here = os.path.dirname(__file__)
    angle_folder = os.path.join(here, '../images/binary/angles/')
    if not os.path.exists(angle_folder):
        os.mkdir(angle_folder)
    try:
        size = int(sys.argv[1])
    except (IndexError, ValueError):
        size = 512

    w = h = size
    b = h / 2  # Height of the triangle
    A = (w / 2, b)  # Tip of the triangle
    for angle in np.linspace(15, 160, 30):
        angle_rad = np.deg2rad(angle / 2)
        name = os.path.join(angle_folder, f'angle{angle:03.0f}-{{}}.pgm')
        B = (w / 2 + b * np.tan(angle_rad), 0)
        C = (w / 2 - b * np.tan(angle_rad), 0)
        im = Image.new('L', (w, h), color=255)
        draw = ImageDraw.Draw(im)
        draw.polygon([C, B, A], fill='black')
        im.save(name.format('0'))
        # im = im.rotate(90)
        # im.save(name.format('1'))
        # im = im.rotate(90)
        # im.save(name.format('2'))
        # im = im.rotate(90)
        # im.save(name.format('3'))
