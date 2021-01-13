"""
This tool generates pseudo H&E images from multidimensional TIFF files containing 
autofluorescence and DAPI channels.

Autofluorescence and DAPI signals used for H&E rendering are maximum intensity projections
of (by default 5 micrometer thick) subsections. 
The location of these subsections are selected for both channels seperately by determining
the z-coordinate of maximum intensity in each channel, but averages over all image tiles.

Current limitations:
- The entire image file (with all tiles) has to fit into RAM.
- 16 bit input images are assumed to contain gray values in the range [0, 2^12-1]


Copyright 2021 Ralf Haeusler

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import sys
import getopt
import numpy as np
from tifffile import imread
from tifffile import imsave
from skimage import exposure


# Color definition for pseudo HnE stain
HnE_COLORS_4DAPI = [0.17, 0.27, 0.105]
HnE_COLORS_4AUTOFLUORESCENCE = [0.05, 1.0, 0.54]


# Histogram equalization potentially useful for application to autofluorescence channel
def histogram_equalize(img):
    img = rgb2gray(img)
    img_cdf, bin_centers = exposure.cumulative_distribution(img)
    return np.interp(img, bin_centers, img_cdf)


# Returns a (by default 5 micron thick) substack, depending on suitable definition of axial image resolution.
def getMinMaxZcoordinatesForSubsection(input5dImage, micronsPerPixZ=0.67, subsectionThicknessInMicron=5):
    # Sum of pixel grayvalues for each z-layer, to obtain
    # maximum position over the z channel. Result for each tile and each channel.
    xySums = np.sum(np.sum(input5dImage, axis=4), axis=3)
    zMaxs = np.argmax(xySums, axis=1)

    # Median z_max over all tiles is center of selected 5um subsection ...
    zMedian = np.median(zMaxs, axis=0)
    zSelected_Dapi = int(zMedian[1])
    print ('Initially selected z coordinate centered at sub-stack: ', zSelected_Dapi)
   
    # ... whose start and end coordinates ...
    zCoordLenHalf = int((subsectionThicknessInMicron/abs(micronsPerPixZ)) / 2)
    layer5um_minZ = zSelected_Dapi - zCoordLenHalf
    layer5um_maxZ = zSelected_Dapi + zCoordLenHalf
    print ('Thickness of selected sub-stack in pixels: ', 2 * zCoordLenHalf + 1)    

    # ... are adjusted if one of the coordinates overextends the z-range of the full stack.
    # In this case, the subsection is shifted as far as possible to fit in 5um.
    maxZ = input5dImage.shape[-4]
    if (layer5um_minZ < 0):
       layer5um_maxZ = layer5um_maxZ - layer5um_minZ
    if (layer5um_maxZ > maxZ):
       layer5um_minZ = layer5um_minZ - (layer5um_maxZ - maxZ)	
     
    # If the ajusted range still overextends the original z-stack, 
    # take the entire z stack, and output warning that layer < 5 um.
    if (layer5um_minZ < 0):
       layer5um_minZ = 0
       print ("Warning: Section selected for maximum intensity projection is less than 5 um thick.")
       
    return layer5um_minZ, layer5um_maxZ
    

def main(argv):
    inputfile = ''
    HNEoutputfile = 'HNEimageRGB_8bit.tif'
    DAPIoutputfile = '' # maxIntensityProjDAPI.tif
    AUTOFLoutputfile = '' # maxIntensityProjAutoFl.tif
    umPerPixZ = 0.67
    subsectionThicknessInMicron = 5
    colorWeight_DAPI = 2.56
    colorWeight_AUTOFLUORESCENCE = 0.1

    try:
       opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print ('Basic usage: python TIFF2HNE.py -i <inputfile> -H <HNE outputfile>')
       print ('Type TIFF2HNE.py -h for additional options.')
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
           print ('python TIFF2HNE.py -i <inputfile> -H <HNE outputfile> ')
           print ('    -A <autofluorescence outputfile> -D <DAPI outputfile>')
           print ('    -z <axial resolution in microns per pixel> (default: ', umPerPixZ, ')')
           print ('    -t <thickness of subsection to select>) (default: ', subsectionThicknessInMicron, ')')
           print ('    -u <color weight for DAPI channel (default: ', colorWeight_DAPI, ')')
           print ('    -v <color weight for autofluorescense channel (default: ', colorWeight_AUTOFLUORESCENCE, ')')
           sys.exit()
       elif opt in ("-i", "--ifile"):
           inputfile = arg
       elif opt in ("-H", "--HNEofile"):
           HNEoutputfile = arg
       elif opt in ("-A", "--AUTOFLofile"):
           AUTOFLoutputfile = arg
       elif opt in ("-D", "--DAPIofile"):
           DAPIoutputfile = arg          
       elif opt in ("-z", "--zRes"):
           umPerPixZ = arg
       elif opt in ("-t", "--subSectionThickness"):          
           subsectionThicknessInMicron = arg
       elif opt in ("-u", "--colWeightDAPI"):
           colorWeight_DAPI = arg
       elif opt in ("-v", "--colWeightAUTOFL"):
           colorWeight_AUTOFLUORESCENCE = arg

    try:
        imageIn = imread(inputfile)
    except FileNotFoundError:
        sys.exit('Input image file does not exist.')                   
            
    if (len(imageIn.shape) < 4 or len(imageIn.shape) > 5):
    	sys.exit("Unsuitable image format: Expected dimension is 4 or 5.")

    if (len(imageIn.shape) == 4):   	   
        print ('Input image dimensions [maxZ, channels, maxY,maxX] and pixel type: ')
    else:
        print ('Input image dimensions [tiles, maxZ, channels, maxY,maxX] and pixel type: ')        
    print(imageIn.shape, ' ',  ' ', imageIn.dtype)    
   
    # Wrap into 5D array if image is only 1 tile.
    imageIn = np.array([imageIn]) if len(imageIn.shape) == 4 else imageIn

    # For each tile cut out 5 um layer. 	  	   
    layer5um_minZ, layer5um_maxZ = getMinMaxZcoordinatesForSubsection(input5dImage=imageIn, micronsPerPixZ=umPerPixZ, \
        subsectionThicknessInMicron=subsectionThicknessInMicron)
    slices5um = imageIn[:,layer5um_minZ:layer5um_maxZ,:,:,:]

    # Maximum intensity projection, contrast normalised, for DAPI and autofluoresence channel.
    # For more contrast adjustment options,
    # see https://scikit-image.org/docs/dev/auto_examples/color_exposure/plot_equalize.html
    maxIntensityProjAUTOFL = np.max(slices5um[:,:,0,:,:], axis=1) 
    maxIntensityProjDAPI = np.max(slices5um[:,:,1,:,:], axis=1) 
    p2AUTOFL, p98AUTOFL = np.percentile(maxIntensityProjAUTOFL, (2, 98))   		   
    p2DAPI, p98DAPI = np.percentile(maxIntensityProjDAPI, (2, 98))   
    contrastAdjustedAUTOFLUORESCENCE = exposure.rescale_intensity(maxIntensityProjAUTOFL, in_range=(p2AUTOFL, p98AUTOFL))    
    contrastAdjustedDAPI = exposure.rescale_intensity(maxIntensityProjDAPI, in_range=(p2DAPI, p98DAPI))
	   	    
    # Compose color image tiles
    RGB_image = np.zeros((contrastAdjustedDAPI.shape[0], 3, contrastAdjustedDAPI.shape[1], contrastAdjustedDAPI.shape[2])) 
    maxGray = 16384 if (slices5um.dtype == np.uint16) else 256
    for tile in range(RGB_image.shape[0]):    
        for rgbChannel in range(RGB_image.shape[1]):
            tmp_AUTOFL = HnE_COLORS_4AUTOFLUORESCENCE[rgbChannel] * colorWeight_AUTOFLUORESCENCE \
                * contrastAdjustedAUTOFLUORESCENCE[tile] / maxGray
            tmp_DAPI = HnE_COLORS_4DAPI[rgbChannel] * colorWeight_DAPI \
                * contrastAdjustedDAPI[tile] / maxGray 
            RGB_image[tile][rgbChannel] = 255 * np.multiply(np.exp(-tmp_AUTOFL), np.exp(-tmp_DAPI))

    # Reshape to layout [tile, Y, X, C] and convert to 24 bit RGB color
    HNE_image = np.moveaxis(RGB_image, 1, -1).astype(np.uint8) 
   
    # For param 'photometric' (below) options in image saving, see: 
    # https://github.com/cgohlke/tifffile/blob/master/tifffile/tifffile.py#L378        
    imsave(HNEoutputfile, HNE_image, photometric='rgb') 
    if DAPIoutputfile != '':
        imsave(DAPIoutputfile, contrastAdjustedDAPI, photometric='minisblack') 
    if AUTOFLoutputfile != '':       
        imsave(AUTOFLoutputfile, contrastAdjustedAUTOFLUORESCENCE, photometric='minisblack')    


if __name__ == "__main__":
   main(sys.argv[1:])
   
