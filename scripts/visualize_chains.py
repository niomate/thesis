from PIL import Image
import re
from matplotlib import colors as cols


def hex_to_rgb(s):
    r = s[1:3]
    g = s[3:5]
    b = s[5:7]
    return (int(r, 16), int(g, 16), int(b, 16))


chain_re = r'([0-9]+), ([0-9]+)'
sort_re = r'angle([0-9][0-9][0-9])-[0-3].pgm'

with open('output.log', 'r') as f:
    lines = f.readlines()

chains = {}
current = ''

for line in lines:
    if '[angle' in line:
        current = line.split(':')[0][1:-1]
        chains[current] = []
    elif 'Current chain' in line:
        matches = re.findall(chain_re, line.split(':')[1])
        chains[current].append([(int(x), int(y)) for (x, y) in matches])

keys = sorted(chains.keys(), key=lambda x: int(re.match(sort_re, x).group(1)))

# TODO: Maybe implement something to go through the different chains and
# visualise them

for key in keys:
    print(key)
    chain = chains[key]
    img = Image.new('RGB', (512, 512), 'white')
    pixels = img.load()
    colors = iter(cols.cnames.values())
    for c in chain:
        current_color = hex_to_rgb(next(colors))
        for x, y in c:
            pixels[x, y] = current_color
    img.show()
    input('Press Enter for next image')
