# coding=utf-8
from PIL import Image
import sys

orig = Image.open(sys.argv[1])
mask = Image.open(sys.argv[2])
outname = sys.argv[2]

if len(sys.argv) > 3:
    outname = sys.argv[3]

size = w, h = orig.size

nmask = Image.new('L', size, color=127)
nmask_pix = nmask.load()

for i in range(w):
    for j in range(h):
        if mask.getpixel((i, j)) == 255:
            nmask_pix[i, j] = orig.getpixel((i, j))

nmask.save(outname)
