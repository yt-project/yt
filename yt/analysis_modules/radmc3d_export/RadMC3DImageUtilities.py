'''

Functions for dealing with the image.out files created by RADMC-3D

'''

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np


def parse_radmc3d_image_header(lines):
    '''

    Parses the header lines from the image.out file.
    Returns a dictionary containing the image parameters.

    '''

    header_data = {}

    # an integer flag describing the image format.
    # this function only works for images made in
    # "observer at infinity" mode, i.e. iformat is 1
    iformat = np.int64(lines[0].strip())
    assert(iformat == 1)
    header_data["iformat"] = iformat

    # The number of pixels in the x and y-directions
    Nx, Ny = [np.int64(Npix) for Npix in lines[1].strip().split()]
    header_data["Nx"] = Nx
    header_data["Ny"] = Ny

    # You can produce images at multiple wavelenths in a single
    # pass. This function assumes that there is only 1.
    num_wavelengths = np.int64(lines[2].strip())
    assert(num_wavelengths == 1)
    header_data["num_wavelengths"] = num_wavelengths

    # The size of pixel in each direction. Note that
    # this only makes sense if iformat is 1
    pixel_size_cm_x, pixel_size_cm_y = \
        [np.float64(Npix) for Npix in lines[3].strip().split()]
    header_data["pixel_size_cm_x"] = pixel_size_cm_x
    header_data["pixel_size_cm_y"] = pixel_size_cm_y

    # The wavelength at which the image was produced.
    # We assume there is only 1 image here.
    wavelength_microns = np.float64(lines[4].strip())  # assume 1 wavelength
    header_data["wavelength_microns"] = wavelength_microns

    return header_data


def read_radmc3d_image(filename):
    '''

    Loads the image.out file created by radmc-3d.
    Returns an np.array that contains the image data
    as well as a dictionary with some useful metadata.

    '''

    fileh = open(filename, 'r')
    lines = fileh.readlines()
    fileh.close()

    # The header should always be 5 lines long,
    # as per the radmc-3d manual
    header_lines = lines[0:6]
    header = parse_radmc3d_image_header(header_lines)

    Nx = header["Nx"]
    Ny = header["Ny"]

    # The rest of the lines are the image data, with the
    # possible exception of a newline at the end
    image_lines = lines[6:]
    image = np.array([np.float64(line.strip()) for line in image_lines
                      if not line.isspace()])

    # This had better be true
    assert(image.size == Nx*Ny)

    image = image.reshape(header["Nx"], header["Ny"])

    return header, image
