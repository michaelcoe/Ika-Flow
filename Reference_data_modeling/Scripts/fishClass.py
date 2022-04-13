from numpy.core.numeric import _cross_dispatcher
import surfaceAreaEstimators as sea
import volumeEstimators as ve

class fish:
    def __init__(self, specimen):
        self.specimen = specimen

        self.sidePolyTop = 0.0
        self.sidePolyBottom = 0.0
        self.nacaFit = 0.0
        self.sideAreaRatio = 0.0
        self.topAreaRatio = 0.0
        self.surfaceArea = 0.0
        self.D_calc = 0.0
        self.volume = 0.0

    def estimate_surface_area(self, topForm, crossSection, aspectRatio, length):
        self.surfaceArea, self.D_calc = sea.determine_surface_area(topForm, crossSection, aspectRatio, length, 
        self.sidePolyTop, self.sidePolyBottom, self.nacaFit, self.nacaFit)

    def estimate_volume(self, topForm, crossSection, aspectRatio, length):
        self.volume = ve.determine_volume(topForm, crossSection, aspectRatio, length,
        self.sidePolyTop, self.sidePolyBottom, self.nacaFit, self.nacaFit)