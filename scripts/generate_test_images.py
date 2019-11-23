from PIL import Image
from PIL import ImageDraw
import random
import os
import sys

NAME_TEMPLATE = os.path.join(os.path.dirname(__file__), '../images/binary/{}.pgm')


def rand_coord():
    return random.randint(1, 255), random.randint(1, 255)


if __name__ == '__main__':
    n_images = int(sys.argv[1])
    for n in range(n_images):
        name = NAME_TEMPLATE.format(f'test{n}')
        im = Image.new('L', (255, 255), color=255)
        points = [rand_coord() for _ in range(4)]
        draw = ImageDraw.Draw(im)
        draw.polygon(points, fill="black")
        im.save(name)

    
