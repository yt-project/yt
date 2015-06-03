cimport cython
cimport numpy as np
import numpy as np
cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_scene as rtcs

def hello_world():
    print "Hello!"

cdef void error_printer(const rtc.RTCError code, const char *_str):
    print "ERROR CAUGHT IN EMBREE"
    rtc.print_error(code)
    print "ERROR MESSAGE:", _str

cdef class EmbreeVolume:

    def __init__(self):
        rtc.rtcInit(NULL)
        rtc.rtcSetErrorFunction(error_printer)
        self.scene_i = rtcs.rtcNewScene(rtcs.RTC_SCENE_STATIC, rtcs.RTC_INTERSECT1)

    def __dealloc__(self):
        rtcs.rtcDeleteScene(self.scene_i)
