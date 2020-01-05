from PIL import Image
from PIL import ImageDraw
import os
import numpy as np

if __name__ == '__main__':
    w, h = 512, 512
    b = 1/2 * h  # Height of the triangle
    A = (w / 2, b)  # Tip of the triangle
    for angle in np.linspace(15, 160, 30):
        angle_rad = np.deg2rad(angle / 2)
        name = os.path.join(
            os.path.dirname(__file__), f'../images/binary/angles/angle{angle:03.0f}-{{}}.pgm')
        B = (w / 2 + b * np.tan(angle_rad), 0)
        C = (w / 2 - b * np.tan(angle_rad), 0)
        im = Image.new('L', (w, h), color=255)
        draw = ImageDraw.Draw(im)
        draw.polygon([C, B, A], fill='black')
        im.save(name.format('0'))
        im = im.rotate(90)
        im.save(name.format('1'))
        im = im.rotate(90)
        im.save(name.format('2'))
        im = im.rotate(90)
        im.save(name.format('3'))
