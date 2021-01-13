# TIFF2HNE.py

Tool to generate pseudo H&E images from multidimensional TIFF files containing 
autofluorescence and DAPI channels.
For 5D images, with the layout `[numTiles, zDim, numChannels, yDim, xDim]`,
the autofluorescence signal is expected to be in `[:,:,0,:,:]`, and the DAPI signal in `[:,:,1,:,:]`.

Autofluorescence and DAPI signals used for H&E rendering are maximum intensity projections
of (by default 5 micrometer thick) subsections. 
The location of these subsections are selected for both channels seperately by determining
the z-coordinate of maximum intensity in each channel, but averages over all image tiles.

Current limitations:
- The entire image file (with all tiles) has to fit into RAM.
- 16 bit input images are assumed to contain gray values in the range `[0, 2^12-1]`.

Usage:
```
python TIFF2HNE.py -i <inputfile> -H <HNE outputfile> 
    -A <autofluorescence outputfile> -D <DAPI outputfile>
    -z <axial resolution in microns per pixel> (default: 5 )
    -t <thickness of subsection to select>) (default: 0.67)
    -u <color weight for DAPI channel (default: 2.56)
    -v <color weight for autofluorescense channel (default: 0.1)
```
