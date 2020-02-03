from PIL import ImageDraw, Image, ImageOps
import os
import sys
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--invert', action='store_true')
parser.add_argument('-r', '--rotate', action='store_true')
parser.add_argument('-t', '--filetype', default='pgm')
parser.add_argument('size', type=int)

if __name__ == '__main__':
    here = os.path.dirname(__file__)
    angle_folder = os.path.join(here, '../images/binary/angles/')

    args = parser.parse_args()

    if not os.path.exists(angle_folder):
        print(f'Directory {angle_folder} not found. Generating a new one...')
        os.mkdir(angle_folder)
    elif len(os.listdir(angle_folder)) > 0:
        print('Deleting previously generated images...')
        for file in os.listdir(angle_folder):
            os.remove(os.path.join(angle_folder, file))

    w = h = args.size
    b = h / 2  # Height of the triangle
    A = (w / 2, b)  # Tip of the triangle
    for angle in np.linspace(15, 160, 30):
        angle_rad = np.deg2rad(angle / 2)
        name = os.path.join(angle_folder, f'angle{angle:03.0f}-{{}}.{args.filetype}')
        B = (w / 2 + b * np.tan(angle_rad), 0)
        C = (w / 2 - b * np.tan(angle_rad), 0)
        im = Image.new('L', (w, h), color=255)
        draw = ImageDraw.Draw(im)
        draw.polygon([C, B, A], fill='black')
        im.save(name.format('0'))
        if args.invert:
            ImageOps.invert(im).save(name.format('0i'))
        elif args.rotate and args.invert:
            im = im.rotate(90)
            im.save(name.format('1'))
            ImageOps.invert(im).save(name.format('1i'))
            im = im.rotate(90)
            im.save(name.format('2'))
            ImageOps.invert(im).save(name.format('2i'))
            im = im.rotate(90)
            im.save(name.format('3'))
            ImageOps.invert(im).save(name.format('3i'))
        elif args.rotate:
            im = im.rotate(90)
            im.save(name.format('1'))
            im = im.rotate(90)
            im.save(name.format('2'))
            im = im.rotate(90)
            im.save(name.format('3'))
