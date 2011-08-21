"""
A light interface to FreeType2

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
cimport numpy as np
cimport cython

cdef int initialized = 0

from stdio cimport fopen, fclose, FILE

cdef extern from "freetype_includes.h":
    # NOTE that size_t might not be int
    ctypedef struct FT_Library:
        pass

    ctypedef struct FT_Bitmap:
        unsigned char *buffer
        int width
        int rows

    ctypedef struct FT_Vector:
        int x
        int y

    ctypedef struct FT_GlyphSlotRec:
        FT_Bitmap bitmap
        int bitmap_left
        int bitmap_top
        FT_Vector advance
    ctypedef FT_GlyphSlotRec* FT_GlyphSlot

    ctypedef struct FT_FaceRec:
        FT_GlyphSlot glyph
    ctypedef FT_FaceRec* FT_Face

    int FT_LOAD_RENDER

    int FT_Init_FreeType( FT_Library *lib)
    int FT_New_Face(FT_Library library, char* filepathname, long face_index,
                    FT_Face *aface)
    int FT_Set_Char_Size(FT_Face face, double width, double height,
                         unsigned dpi_w, unsigned dpi_h)
    int FT_Load_Char(FT_Face face, unsigned long char_code, int load_flags)

cdef FT_Library library

def initialize_library():
    # We do not clear the library from memory.
    global initialized
    if initialized == 1: return
    cdef int error = FT_Init_FreeType(&library)
    if error: raise RuntimeError
    initialized = 1

def simple_writing(char *fontpath, long face_index, long dpi,
                   long font_size, char *text,
                   np.ndarray[np.uint8_t, ndim=3] buffer,
                   int xpos, int ypos):
    initialize_library()
    cdef FT_Face face
    cdef int error = FT_New_Face(library, fontpath, face_index, &face)
    if error: raise RuntimeError(error)
    error = FT_Set_Char_Size(face, 0, font_size * 64, dpi, dpi)
    if error: raise RuntimeError(error)
    cdef FT_GlyphSlot slot = face.glyph
    cdef int n, i, j, c
    cdef pi = xpos, pj = ypos # Pen locations
    cdef di, dj # Where we actually draw
    cdef char tchar
    cdef int ichar
    cdef np.uint8_t bb
    for n in range(len(text)):
        tchar = text[n]
        ichar = <int> tchar
        error = FT_Load_Char(face, ichar, FT_LOAD_RENDER)
        if error: raise RuntimeError(error)
        for i in range(slot.bitmap.width):
            for j in range(slot.bitmap.rows):
                di = pi + slot.bitmap_left + i
                dj = pj - slot.bitmap_top + j
                for c in range(3):
                    bb = <np.uint8_t> slot.bitmap.buffer[j * slot.bitmap.width + i]
                    buffer[dj,di,c] |= bb
        pi += slot.advance.x / 64
        pj += slot.advance.y / 64
