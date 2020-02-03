import subprocess
import os

path = '../images/binary/angles/'
for file in os.listdir(path):
    cmd = ['bin/harris_corner_detector', os.path.join(path, file), '-v', '-o', f'{file.split(".")[0]}_corners.png']
    print(f'Running {" ".join(cmd)}')
    subprocess.run(cmd)

