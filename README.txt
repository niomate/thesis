######################################################################
#                                                                    #
#         README - Usage for corner_detection and inpaint            #
#                                                                    #
#                  Daniel Gusenburger July 2020                      #
#                                                                    #
######################################################################

To build both programs, simply run make.

1. Corner detection with the structure tensor

./corner [OPTIONS] [input image file (.pgm format)]

Options:
--------

	-s <float>:	Sigma/Noise scale 
			Default: 1.0

	-r <float>:	Rho/Integration scale
			Default: 2.5

	-c <int>:     	Type of the corner detector
			  0 -> Rohr
			  1 -> Tomasi-Kanade
			  2 -> Förstner-Harris
			Default: Förstner-Harris (2)
	
	-q <float>:	Percentile parameter 
			Meaning depends on whether TPPT is enabled or not
			Default: 0.1

	-m <float>:     Mask radius in pixels
			Default: 0

	-C:		Flag to enable CNMS, uses -m as mask radius

	-P:		Flag to enable TPPT, uses -q as percentile

	-M: 		Output the computed mask (in black/white binary) instead of the corner location

	-d <char>:      'x' to mark corner locations as crosses
			'o' to mark corner locations as circles
			Default: 'o'

	-o <string>:	Filename for output file.
			Default: "corners.pgm"


2. Inpainting using edge-enhancing diffusion

######################################################################
#                       START DISCLAIMER                             # 
######################################################################
#                                                                    #
# The mask image used in this version has to be a black white image, #
# where the seed points are marked white and the rest is black.      #
# The other input image has to be the original image.                #
# Using this mask and the original image, the initial image is       #
# created as an image where the mask points are initiated by their   #
# value in the original image and all other values are initiated     #
# with a grey value of 127.                                          #
#                                                                    #
######################################################################
#                       END DISCLAIMER                               # 
######################################################################

./inpaint [OPTIONS] [input image - original (.pgm format)] [mask image (.pgm format)]

Options:
--------

	-s <float>:	Sigma/Noise scale
			Default: 3.0

	-l <float>: 	Lambda/Contrast parameter
			Default: 1.0

	-g <float>: 	Gamma/Non-negativity parameter
			Default: 1.0

	-a <float>:	Alpha/Dissipitavity parameter
			Default: 0.49
	
	-n <int>:	Maximum number of outer iterations
			Default: 1000

	-N <int>:       Number of inner/solver iterations
			Default: 200

	-o <string>:	Filename of output image
			Default: "inpaint.pgm"

As for the other parameters in the original version:

Diffusivity: 		Charbonnier
Time-discretisation:	semi-implicit
Solver:			Conjugate Gradients (CG)

Rho is always set to 0 since I examined EED not CED.
