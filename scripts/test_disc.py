from run_all import corners, inpaint, mask, name_from_path
import pathlib

params_inpaint = {'-s': 3.0, '-l': 1.0, '-a': 0.49,
                  '-g': 1.0, '-n': '1000', '-N': '200'}
params_same = {'-q': 0.01, '-s': 1, '-r': 4, '-m': 1, '-c': 1, '-D': ''}
params_mult1 = {'-q': 0.01, '-s': 1, '-r': 2, '-m': 2, '-c': 1, '-D': ''}
params_mult2 = {'-q': 0.01, '-s': 1, '-r': 1, '-m': 4, '-c': 1, '-D': ''}

params = [params_same, params_mult1, params_mult2]


def flatten(dic):
    return [x for kv in dic.items() for x in kv]


path = pathlib.Path('~/Uni/Thesis/').expanduser()

image_dir = path / 'images' / 'binary' / 'angles'
output_dir = path / 'output' / 'disctest'

if not output_dir.exists():
    output_dir.mkdir()

for img in image_dir.iterdir():
    n, _ = img.name.rsplit('.', 1)
    for p in params:
        out_name_corners = f'{n}{p["-r"]}{p["-m"]}corners.pgm'
        out_name_mask = f'{n}{p["-r"]}{p["-m"]}mask.pgm'
        out_name_inpaint = f'{n}{p["-r"]}{p["-m"]}inpaint.pgm'

        outim = output_dir / out_name_corners
        corners(img.as_posix(), *flatten(p), outim=outim.as_posix())
        outim_mask = output_dir / out_name_mask
        mask(img.as_posix(), *flatten(p), outim=outim_mask.as_posix())
        outim = output_dir / out_name_inpaint
        inpaint(img.as_posix(), outim_mask.as_posix(), *
                flatten(params_inpaint), outim=outim.as_posix())
