######################################### PYTHON SCRIPTS #################################################

Usage for run.py:

python run.py [ corners | mask | full ] [corners command line options]

run.py essentially runs the corner method for all .pgm images found in images/binary concurrently.

It can be given all command line arguments specified below in the corners usage section.
Output images are stored in output/masks and output/inpaint (only if 'full' was specified)

######################################### CORNER DETECTION ###############################################

Usage for corner detection: 

corners [OPTIONS] input_image

OPTIONS:

    -q: percentile ratio for corner detection (0 <= q <= 1)

    -c: measure type
        0: Rohr
        1: Tomasi-Kanade
        2: Foerstner-Harris
        3: Alternate Harris measure
    
    -s sigma for preprocessing gaussian smoothing

    -r sigma for integration scale

    -k kappa for alternate harris measure

    -o output output image (if not specified, saves the result in "corners.pgm")

    -m mask radius
        If specified, creates a mask image for corner inpainting with the given radius


############################################ INPAINTING ##################################################

Usage for inpainting:

Soma parameters are assumed and can't be changed via command line, 
such as time discretization and the solver method.


inpainting [OPTIONS] input_image mask

OPTIONS:

    -l lambda 

    -s sigma

    -a alpha

    -g gamma
    
    -n outer iterations 

    -N solver iterations

    -o output image (if not specified, saves the result in "inpaint.pgm")
