"""


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np

import matplotlib.colors as cc
import matplotlib.cm as mcm
from . import _colormap_data as _cm
from yt.extern.six import string_types

try:
    import cmocean
except ImportError:
    cmocean = None

def is_colormap(cmap):
    return isinstance(cmap,cc.Colormap)

def check_color(name):
    try:
        cc.colorConverter.to_rgb(name)
        return True
    except ValueError:
        return False

yt_colormaps = {}

def add_cmap(name, cdict):
    """Deprecated alias, kept for backwards compatibility."""
    from yt.funcs import issue_deprecation_warning
    issue_deprecation_warning("Deprecated alias. Use add_colormap instead.")
    add_colormap(name, cdict)

def add_colormap(name, cdict):
    """
    Adds a colormap to the colormaps available in yt for this session
    """
    yt_colormaps[name] = \
        cc.LinearSegmentedColormap(name,cdict,256)
    mcm.datad[name] = cdict
    mcm.__dict__[name] = cdict
    try: # API compatibility
        mcm.register_cmap(name, yt_colormaps[name])
    except AttributeError:
        pass
    

# The format is as follows:
#   First number is the number at which we are defining a color breakpoint
#   Second number is the (0..1) number to interpolate to when coming *from below*
#   Third number is the (0..1) number to interpolate to when coming *from above*

# Next up is boilerplate -- the name, the colormap dict we just made, and the
# number of segments we want.  This is probably fine as is.

cdict = {'red':   ((0.0, 80/256., 80/256.),
                   (0.2, 0.0, 0.0),
                   (0.4, 0.0, 0.0),
                   (0.6, 256/256., 256/256.),
                   (0.95, 256/256., 256/256.),
                   (1.0, 150/256., 150/256.)),
         'green': ((0.0, 0/256., 0/256.),
                   (0.2, 0/256., 0/256.),
                   (0.4, 130/256., 130/256.),
                   (0.6, 256/256., 256/256.),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 80/256., 80/256.),
                   (0.2, 220/256., 220/256.),
                   (0.4, 0.0, 0.0),
                   (0.6, 20/256., 20/256.),
                   (1.0, 0.0, 0.0))}

add_colormap('bds_highcontrast', cdict)
add_colormap('algae', cdict)

# This next colormap was designed by Tune Kamae and converted here by Matt
_vs = np.linspace(0,1,255)
_kamae_red = np.minimum(255,
                113.9*np.sin(7.64*(_vs**1.705)+0.701)-916.1*(_vs+1.755)**1.862 \
              + 3587.9*_vs+2563.4)/255.0
_kamae_grn = np.minimum(255,
                70.0*np.sin(8.7*(_vs**1.26)-2.418)+151.7*_vs**0.5+70.0)/255.0
_kamae_blu = np.minimum(255,
                194.5*_vs**2.88+99.72*np.exp(-77.24*(_vs-0.742)**2.0)
              + 45.40*_vs**0.089+10.0)/255.0

cdict = {'red':np.transpose([_vs,_kamae_red,_kamae_red]),
         'green':np.transpose([_vs,_kamae_grn,_kamae_grn]),
         'blue':np.transpose([_vs,_kamae_blu,_kamae_blu])}
add_colormap('kamae', cdict)

# This one is a simple black & green map

cdict = {'red':   ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),
         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0))}

add_colormap('black_green', cdict)

cdict = {'red':   ((0.0, 0.0, 0.0),
                   (1.0, 0.2, 0.2)),
         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.2, 0.2)),
         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0))}

add_colormap('black_blueish', cdict)

# This one is a variant of a colormap commonly
# used for X-ray observations by Maxim Markevitch

cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.3, 0.0, 0.0),
                 (0.352, 0.245, 0.245),
                 (0.42, 0.5, 0.5),
                 (0.51, 0.706, 0.706),
                 (0.613, 0.882, 0.882),
                 (0.742, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),
                   (0.585, 0.0, 0.0),
                   (0.613, 0.196, 0.196),
                   (0.693, 0.48, 0.48),
                   (0.785, 0.696, 0.696),
                   (0.885, 0.882, 0.882),
                   (1.0, 1.0, 1.0)),
         'blue': ((0.0, 0.0, 0.0),
                  (0.136, 0.0, 0.0),
                  (0.136, 0.373, 0.373),
                  (0.391, 1.0, 1.0),
                  (1.0, 1.0, 1.0))}

add_colormap("purple_mm", cdict)

# This one comes from
# http://permalink.gmane.org/gmane.comp.python.matplotlib.devel/10518
# and is an implementation of http://arxiv.org/abs/1108.5083
#

# cubehelix parameters
_gamma_cubehelix = 1.0
_s_cubehelix = 0.5
_r_cubehelix = -1.5
_h_cubehelix = 1.0

_cubehelix_data = {
        'red': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (-0.14861 * np.cos(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x)) + 1.78277 * np.sin(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
        'green': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (-0.29227 * np.cos(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x)) - 0.90649 * np.sin(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
        'blue': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (1.97294 * np.cos(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
}

add_colormap("cubehelix", _cubehelix_data)

# The turbo colormap, by Anton Mikhailov.
# from: https://gist.github.com/mikhailov-work/ee72ba4191942acecc03fe6da94fc73f

_turbo_colormap_data = np.array(
    [[0.18995,0.07176,0.23217],[0.19483,0.08339,0.26149],
     [0.19956,0.09498,0.29024],[0.20415,0.10652,0.31844],
     [0.20860,0.11802,0.34607],[0.21291,0.12947,0.37314],
     [0.21708,0.14087,0.39964],[0.22111,0.15223,0.42558],
     [0.22500,0.16354,0.45096],[0.22875,0.17481,0.47578],
     [0.23236,0.18603,0.50004],[0.23582,0.19720,0.52373],
     [0.23915,0.20833,0.54686],[0.24234,0.21941,0.56942],
     [0.24539,0.23044,0.59142],[0.24830,0.24143,0.61286],
     [0.25107,0.25237,0.63374],[0.25369,0.26327,0.65406],
     [0.25618,0.27412,0.67381],[0.25853,0.28492,0.69300],
     [0.26074,0.29568,0.71162],[0.26280,0.30639,0.72968],
     [0.26473,0.31706,0.74718],[0.26652,0.32768,0.76412],
     [0.26816,0.33825,0.78050],[0.26967,0.34878,0.79631],
     [0.27103,0.35926,0.81156],[0.27226,0.36970,0.82624],
     [0.27334,0.38008,0.84037],[0.27429,0.39043,0.85393],
     [0.27509,0.40072,0.86692],[0.27576,0.41097,0.87936],
     [0.27628,0.42118,0.89123],[0.27667,0.43134,0.90254],
     [0.27691,0.44145,0.91328],[0.27701,0.45152,0.92347],
     [0.27698,0.46153,0.93309],[0.27680,0.47151,0.94214],
     [0.27648,0.48144,0.95064],[0.27603,0.49132,0.95857],
     [0.27543,0.50115,0.96594],[0.27469,0.51094,0.97275],
     [0.27381,0.52069,0.97899],[0.27273,0.53040,0.98461],
     [0.27106,0.54015,0.98930],[0.26878,0.54995,0.99303],
     [0.26592,0.55979,0.99583],[0.26252,0.56967,0.99773],
     [0.25862,0.57958,0.99876],[0.25425,0.58950,0.99896],
     [0.24946,0.59943,0.99835],[0.24427,0.60937,0.99697],
     [0.23874,0.61931,0.99485],[0.23288,0.62923,0.99202],
     [0.22676,0.63913,0.98851],[0.22039,0.64901,0.98436],
     [0.21382,0.65886,0.97959],[0.20708,0.66866,0.97423],
     [0.20021,0.67842,0.96833],[0.19326,0.68812,0.96190],
     [0.18625,0.69775,0.95498],[0.17923,0.70732,0.94761],
     [0.17223,0.71680,0.93981],[0.16529,0.72620,0.93161],
     [0.15844,0.73551,0.92305],[0.15173,0.74472,0.91416],
     [0.14519,0.75381,0.90496],[0.13886,0.76279,0.89550],
     [0.13278,0.77165,0.88580],[0.12698,0.78037,0.87590],
     [0.12151,0.78896,0.86581],[0.11639,0.79740,0.85559],
     [0.11167,0.80569,0.84525],[0.10738,0.81381,0.83484],
     [0.10357,0.82177,0.82437],[0.10026,0.82955,0.81389],
     [0.09750,0.83714,0.80342],[0.09532,0.84455,0.79299],
     [0.09377,0.85175,0.78264],[0.09287,0.85875,0.77240],
     [0.09267,0.86554,0.76230],[0.09320,0.87211,0.75237],
     [0.09451,0.87844,0.74265],[0.09662,0.88454,0.73316],
     [0.09958,0.89040,0.72393],[0.10342,0.89600,0.71500],
     [0.10815,0.90142,0.70599],[0.11374,0.90673,0.69651],
     [0.12014,0.91193,0.68660],[0.12733,0.91701,0.67627],
     [0.13526,0.92197,0.66556],[0.14391,0.92680,0.65448],
     [0.15323,0.93151,0.64308],[0.16319,0.93609,0.63137],
     [0.17377,0.94053,0.61938],[0.18491,0.94484,0.60713],
     [0.19659,0.94901,0.59466],[0.20877,0.95304,0.58199],
     [0.22142,0.95692,0.56914],[0.23449,0.96065,0.55614],
     [0.24797,0.96423,0.54303],[0.26180,0.96765,0.52981],
     [0.27597,0.97092,0.51653],[0.29042,0.97403,0.50321],
     [0.30513,0.97697,0.48987],[0.32006,0.97974,0.47654],
     [0.33517,0.98234,0.46325],[0.35043,0.98477,0.45002],
     [0.36581,0.98702,0.43688],[0.38127,0.98909,0.42386],
     [0.39678,0.99098,0.41098],[0.41229,0.99268,0.39826],
     [0.42778,0.99419,0.38575],[0.44321,0.99551,0.37345],
     [0.45854,0.99663,0.36140],[0.47375,0.99755,0.34963],
     [0.48879,0.99828,0.33816],[0.50362,0.99879,0.32701],
     [0.51822,0.99910,0.31622],[0.53255,0.99919,0.30581],
     [0.54658,0.99907,0.29581],[0.56026,0.99873,0.28623],
     [0.57357,0.99817,0.27712],[0.58646,0.99739,0.26849],
     [0.59891,0.99638,0.26038],[0.61088,0.99514,0.25280],
     [0.62233,0.99366,0.24579],[0.63323,0.99195,0.23937],
     [0.64362,0.98999,0.23356],[0.65394,0.98775,0.22835],
     [0.66428,0.98524,0.22370],[0.67462,0.98246,0.21960],
     [0.68494,0.97941,0.21602],[0.69525,0.97610,0.21294],
     [0.70553,0.97255,0.21032],[0.71577,0.96875,0.20815],
     [0.72596,0.96470,0.20640],[0.73610,0.96043,0.20504],
     [0.74617,0.95593,0.20406],[0.75617,0.95121,0.20343],
     [0.76608,0.94627,0.20311],[0.77591,0.94113,0.20310],
     [0.78563,0.93579,0.20336],[0.79524,0.93025,0.20386],
     [0.80473,0.92452,0.20459],[0.81410,0.91861,0.20552],
     [0.82333,0.91253,0.20663],[0.83241,0.90627,0.20788],
     [0.84133,0.89986,0.20926],[0.85010,0.89328,0.21074],
     [0.85868,0.88655,0.21230],[0.86709,0.87968,0.21391],
     [0.87530,0.87267,0.21555],[0.88331,0.86553,0.21719],
     [0.89112,0.85826,0.21880],[0.89870,0.85087,0.22038],
     [0.90605,0.84337,0.22188],[0.91317,0.83576,0.22328],
     [0.92004,0.82806,0.22456],[0.92666,0.82025,0.22570],
     [0.93301,0.81236,0.22667],[0.93909,0.80439,0.22744],
     [0.94489,0.79634,0.22800],[0.95039,0.78823,0.22831],
     [0.95560,0.78005,0.22836],[0.96049,0.77181,0.22811],
     [0.96507,0.76352,0.22754],[0.96931,0.75519,0.22663],
     [0.97323,0.74682,0.22536],[0.97679,0.73842,0.22369],
     [0.98000,0.73000,0.22161],[0.98289,0.72140,0.21918],
     [0.98549,0.71250,0.21650],[0.98781,0.70330,0.21358],
     [0.98986,0.69382,0.21043],[0.99163,0.68408,0.20706],
     [0.99314,0.67408,0.20348],[0.99438,0.66386,0.19971],
     [0.99535,0.65341,0.19577],[0.99607,0.64277,0.19165],
     [0.99654,0.63193,0.18738],[0.99675,0.62093,0.18297],
     [0.99672,0.60977,0.17842],[0.99644,0.59846,0.17376],
     [0.99593,0.58703,0.16899],[0.99517,0.57549,0.16412],
     [0.99419,0.56386,0.15918],[0.99297,0.55214,0.15417],
     [0.99153,0.54036,0.14910],[0.98987,0.52854,0.14398],
     [0.98799,0.51667,0.13883],[0.98590,0.50479,0.13367],
     [0.98360,0.49291,0.12849],[0.98108,0.48104,0.12332],
     [0.97837,0.46920,0.11817],[0.97545,0.45740,0.11305],
     [0.97234,0.44565,0.10797],[0.96904,0.43399,0.10294],
     [0.96555,0.42241,0.09798],[0.96187,0.41093,0.09310],
     [0.95801,0.39958,0.08831],[0.95398,0.38836,0.08362],
     [0.94977,0.37729,0.07905],[0.94538,0.36638,0.07461],
     [0.94084,0.35566,0.07031],[0.93612,0.34513,0.06616],
     [0.93125,0.33482,0.06218],[0.92623,0.32473,0.05837],
     [0.92105,0.31489,0.05475],[0.91572,0.30530,0.05134],
     [0.91024,0.29599,0.04814],[0.90463,0.28696,0.04516],
     [0.89888,0.27824,0.04243],[0.89298,0.26981,0.03993],
     [0.88691,0.26152,0.03753],[0.88066,0.25334,0.03521],
     [0.87422,0.24526,0.03297],[0.86760,0.23730,0.03082],
     [0.86079,0.22945,0.02875],[0.85380,0.22170,0.02677],
     [0.84662,0.21407,0.02487],[0.83926,0.20654,0.02305],
     [0.83172,0.19912,0.02131],[0.82399,0.19182,0.01966],
     [0.81608,0.18462,0.01809],[0.80799,0.17753,0.01660],
     [0.79971,0.17055,0.01520],[0.79125,0.16368,0.01387],
     [0.78260,0.15693,0.01264],[0.77377,0.15028,0.01148],
     [0.76476,0.14374,0.01041],[0.75556,0.13731,0.00942],
     [0.74617,0.13098,0.00851],[0.73661,0.12477,0.00769],
     [0.72686,0.11867,0.00695],[0.71692,0.11268,0.00629],
     [0.70680,0.10680,0.00571],[0.69650,0.10102,0.00522],
     [0.68602,0.09536,0.00481],[0.67535,0.08980,0.00449],
     [0.66449,0.08436,0.00424],[0.65345,0.07902,0.00408],
     [0.64223,0.07380,0.00401],[0.63082,0.06868,0.00401],
     [0.61923,0.06367,0.00410],[0.60746,0.05878,0.00427],
     [0.59550,0.05399,0.00453],[0.58336,0.04931,0.00486],
     [0.57103,0.04474,0.00529],[0.55852,0.04028,0.00579],
     [0.54583,0.03593,0.00638],[0.53295,0.03169,0.00705],
     [0.51989,0.02756,0.00780],[0.50664,0.02354,0.00863],
     [0.49321,0.01963,0.00955],[0.47960,0.01583,0.01055]])

_tvals = np.linspace(0, 1, 256)
_turbo_data = \
  dict((color, np.transpose([_tvals,
                             _turbo_colormap_data[:, i],
                             _turbo_colormap_data[:, i]]))
       for i, color in enumerate(['red', 'green', 'blue']))

add_colormap("turbo", _turbo_data)

# Add colormaps from cmocean, if it's installed
if cmocean is not None:
    cmo_cmapnames = cmocean.cm.cmapnames
    cmo_cmapnames += ["%s_r" % name for name in cmo_cmapnames]
    for cmname in cmo_cmapnames:
        cm = getattr(cmocean.cm, cmname)
        # cmocean has a colormap named 'algae', so let's avoid overwriting
        # yt's algae or any other colormap we've already added
        if cmname in yt_colormaps:
            cmname = cmname + '_cmocean'
        yt_colormaps[cmname] = cm
        try:
            mcm.register_cmap(cmname, yt_colormaps[cmname])
        except AttributeError:
            # for old versions of matplotlib this won't work, so we avoid
            # erroring out but don't worry about registering with matplotlib
            pass

# Add colormaps in _colormap_data.py that weren't defined here
_vs = np.linspace(0,1,256)
for k,v in list(_cm.color_map_luts.items()):
    if k not in yt_colormaps and k not in mcm.cmap_d:
        cdict = { 'red': np.transpose([_vs,v[0],v[0]]),
                  'green': np.transpose([_vs,v[1],v[1]]),
                  'blue': np.transpose([_vs,v[2],v[2]]) }
        add_colormap(k, cdict)

def _extract_lookup_table(cmap_name):
    cmap = mcm.get_cmap(cmap_name)
    if not cmap._isinit: cmap._init()
    r = cmap._lut[:-3, 0]
    g = cmap._lut[:-3, 1]
    b = cmap._lut[:-3, 2]
    a = np.ones(b.shape)
    return [r, g, b, a]

def show_colormaps(subset="all", filename=None):
    """
    Displays the colormaps available to yt.  Note, most functions can use
    both the matplotlib and the native yt colormaps; however, there are 
    some special functions existing within image_writer.py (e.g. write_image()
    write_bitmap(), etc.), which cannot access the matplotlib
    colormaps.

    In addition to the colormaps listed, one can access the reverse of each 
    colormap by appending a "_r" to any map.

    If you wish to only see certain colormaps, include them in the cmap_list
    attribute.
    
    Parameters
    ----------

    subset : string, or list of strings, optional

        valid values : "all", "yt_native", or list of cmap names
        default : "all"

        As mentioned above, a few functions can only access yt_native 
        colormaps.  To display only the yt_native colormaps, set this
        to "yt_native".  

        If you wish to only see a few colormaps side by side, you can
        include them as a list of colormap names.  
        Example: ['algae', 'gist_stern', 'kamae', 'spectral']

    filename : string, opt

        default: None

        If filename is set, then it will save the colormaps to an output
        file.  If it is not set, it will "show" the result interactively.
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    a=np.outer(np.arange(0,1,0.01), np.ones(10))
    if subset == "all":
        maps = [ m for m in cm.cmap_d if (not m.startswith("idl")) & (not m.endswith("_r"))]
    elif subset == "yt_native":
        maps = [ m for m in _cm.color_map_luts if (not m.startswith("idl")) & (not m.endswith("_r"))]
    else:
        try:
            maps = [ m for m in cm.cmap_d if m in subset]
            if len(maps) == 0:
                raise AttributeError
        except AttributeError:
            raise AttributeError("show_colormaps requires subset attribute "
                                 "to be 'all', 'yt_native', or a list of "
                                 "valid colormap names.")
    maps = list(set(maps))
    maps.sort()
    # scale the image size by the number of cmaps
    plt.figure(figsize=(2.*len(maps)/10.,6))
    plt.subplots_adjust(top=0.7,bottom=0.05,left=0.01,right=0.99)
    l = len(maps)+1
    for i,m in enumerate(maps):
        plt.subplot(1,l,i+1)
        plt.axis("off")
        plt.imshow(a, aspect='auto',cmap=plt.get_cmap(m),origin="lower")      
        plt.title(m,rotation=90, fontsize=10, verticalalignment='bottom')
    if filename is not None:
        plt.savefig(filename, dpi=100, facecolor='gray') 
    else:  
        plt.show()

def make_colormap(ctuple_list, name=None, interpolate=True):
    """
    This generates a custom colormap based on the colors and spacings you
    provide.  Enter a ctuple_list, which consists of tuples of (color, spacing)
    to return a colormap appropriate for use in yt.  If you specify a 
    name, it will automatically be added to the current session as a valid
    colormap.
    
    Output colormap is in the format yt expects for adding a colormap to the 
    current session: a dictionary with the appropriate RGB channels each 
    consisting of a 256x3 array :
    First number is the number at which we are defining a color breakpoint
    Second number is the (0..1) number to interpolate to when coming *from below*
    Third number is the (0..1) number to interpolate to when coming *from above*

    Parameters
    ----------

    ctuple_list: list of (color, float) tuples
        The ctuple_list consists of pairs of (color, interval) tuples
        identifying the colors to use in the colormap and the intervals
        they take to change to the next color in the list.  A color can
        either be a string of the name of a color, or it can be an array
        of 3 floats, each representing the intensity of R, G, and B on
        a scale of 0 to 1.  Valid color names and their equivalent
        arrays are listed below.

        Any interval can be given for the different color tuples, and
        the total of all the intervals will be scaled to the 256 output
        elements.

        If a ctuple_list ends with a color and a non-zero interval, 
        a white 0-interval would be added to the end to finish the 
        interpolation.  To avoid finishing with white, specify your own 
        zero-interval color at the end.

    name: string, optional
        If you wish this colormap to be added as a valid colormap to the 
        current session, specify a name here.  Default: None

    interpolation: boolean, optional
        Designates whether or not the colormap will interpolate between
        the colors provided or just give solid colors across the intervals.
        Default: True

    Preset Color Options
    --------------------

    'white' : np.array([255, 255, 255 ])/255.
    'gray' : np.array([130, 130, 130])/255.
    'dgray' : np.array([80, 80, 80])/255.
    'black' : np.array([0, 0, 0])/255.
    'blue' : np.array([0, 0, 255])/255.
    'dblue' : np.array([0, 0, 160])/255.
    'purple' : np.array([100, 0, 200])/255.
    'dpurple' : np.array([66, 0, 133])/255.
    'dred' : np.array([160, 0, 0])/255.
    'red' : np.array([255, 0, 0])/255.
    'orange' : np.array([255, 128, 0])/255.
    'dorange' : np.array([200,100, 0])/255.
    'yellow' : np.array([255, 255, 0])/255.
    'dyellow' : np.array([200, 200, 0])/255.
    'green' : np.array([0, 255, 0])/255.
    'dgreen' : np.array([0, 160, 0])/255.

    Examples
    --------

    To obtain a colormap that starts at black with equal intervals in green, 
    blue, red, yellow in that order and interpolation between those colors.  
    (In reality, it starts at black, takes an interval of 10 to interpolate to 
    green, then an interval of 10 to interpolate to blue, then an interval of 
    10 to interpolate to red.)

    >>> cm = make_colormap([('black', 10), ('green', 10), ('blue', 10),
    ...                     ('red', 0)])

    To add a colormap that has five equal blocks of solid major colors to
    the current session as "steps":

    >>> make_colormap([('red', 10), ('orange', 10), ('yellow', 10),
    ...                ('green', 10), ('blue', 10)], name="steps",
    ...               interpolate=False)

    To add a colormap that looks like the French flag (i.e. equal bands of 
    blue, white, and red) using your own RGB keys, then to display it:

    >>> make_colormap([([0,0,1], 10), ([1,1,1], 10), ([1,0,0], 10)], 
    ...               name='french_flag', interpolate=False)
    >>> show_colormaps(['french_flag'])

    """
    # aliases for different colors
    color_dict = {
    'white' : np.array([255, 255, 255 ])/255.,
    'gray' : np.array([130, 130, 130])/255.,
    'dgray' : np.array([80, 80, 80])/255.,
    'black' : np.array([0, 0, 0])/255.,
    'blue' : np.array([0, 0, 255])/255.,
    'dblue' : np.array([0, 0, 160])/255.,
    'purple' : np.array([100, 0, 200])/255.,
    'dpurple' : np.array([66, 0, 133])/255.,
    'dred' : np.array([160, 0, 0])/255.,
    'red' : np.array([255, 0, 0])/255.,
    'orange' : np.array([255, 128, 0])/255.,
    'dorange' : np.array([200,100, 0])/255.,
    'yellow' : np.array([255, 255, 0])/255.,
    'dyellow' : np.array([200, 200, 0])/255.,
    'green' : np.array([0, 255, 0])/255.,
    'dgreen' : np.array([0, 160, 0])/255.}

    cmap = np.zeros((256,3))

    # If the user provides a list with a non-zero final interval, it
    # doesn't make sense because you have an interval but no final
    # color to which it interpolates.  So provide a 0-length white final
    # interval to end the previous interval in white.
    if ctuple_list[-1][1] != 0:
        ctuple_list.append(('white', 0))

    # Figure out how many intervals there are total.
    rolling_index = 0
    for i, (color, interval) in enumerate(ctuple_list):
        if isinstance(color, string_types):
            ctuple_list[i] = (color_dict[color], interval)
        rolling_index += interval
    scale = 256./rolling_index
    n = len(ctuple_list)

    # Step through each ctuple and interpolate from one color to the
    # next over the interval provided
    rolling_index = 0
    for i in range(n-1):
        color, interval = ctuple_list[i]
        interval *= scale
        next_index = rolling_index + interval
        next_color, next_interval = ctuple_list[i+1]

        if not interpolate:
            next_color = color

        # Interpolate the R, G, and B channels from one color to the next
        # Use np.round to make sure you're on a discrete index
        interval = int(np.round(next_index)-np.round(rolling_index))
        for j in np.arange(3):
            cmap[int(np.rint(rolling_index)):int(np.rint(next_index)), j] = \
                np.linspace(color[j], next_color[j], num=interval)

        rolling_index = next_index

    # Return a dictionary with the appropriate RGB channels each consisting of
    # a 256x3 array in the format that is expected by add_colormap() to add a 
    # colormap to the session.

    # The format is as follows:
    #   First number is the number at which we are defining a color breakpoint
    #   Second number is the (0..1) number to interpolate to when coming *from below*
    #   Third number is the (0..1) number to interpolate to when coming *from above*
    _vs = np.linspace(0,1,256)
    cdict = {'red':   np.transpose([_vs, cmap[:,0], cmap[:,0]]),
             'green': np.transpose([_vs, cmap[:,1], cmap[:,1]]),
             'blue':  np.transpose([_vs, cmap[:,2], cmap[:,2]])}

    if name is not None:
        add_colormap(name, cdict)

    return cdict
