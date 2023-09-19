#----------------------------------------------------------------------------------------#
#---------------------------------------//GEESEBAL//-------------------------------------#
#GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
#CREATE BY: LEONARDO LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF AND AYAN FLEISCHMANN
#PROJECT - ET BRASIL https://etbrasil.org/
#LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
#UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
#RIO GRANDE DO SUL, BRAZIL

#DOI
#VERSION 0.1.1
#CONTACT US: leonardo.laipelt@ufrgs.br

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#PYTHON PACKAGES
#Call EE
import ee
#ee.Initialize()

#FOLDERS
from .masks import (
f_cloudMaskL457_SR,f_cloudMaskL8_SR,
 f_albedoL5L7,f_albedoL8)
from .meteorology import get_meteorology
from .tools import (fexp_spec_ind, fexp_lst_export,fexp_radlong_up, LST_DEM_correction,
fexp_radshort_down, fexp_radlong_down, fexp_radbalance, fexp_soil_heat,fexp_sensible_heat_flux)
from .endmembers import fexp_cold_pixel, fexp_hot_pixel
from .evapotranspiration import fexp_et

#SEBAL algorithm function
def sebal(NDVI_cold=5, Ts_cold=20, NDVI_hot=10, Ts_hot=20):
    """SEBAL algorithm to be mapped to a collection.
    """
    def fsebal(image):
        p_top_NDVI=ee.Number(NDVI_cold) 
        p_coldest_Ts=ee.Number(Ts_cold)
        p_lowest_NDVI=ee.Number(NDVI_hot)
        p_hottest_Ts=ee.Number(Ts_hot)

        geometryReducer=image.geometry().bounds()
        sun_elevation=ee.Number(90).subtract(image.get("SOLAR_ZENITH_ANGLE"))
        col_meteorology=get_meteorology(ee.ImageCollection(image),image.get("system:time_start"))
        T_air = col_meteorology.select("AirT_G") 
        ux = col_meteorology.select("ux_G")      
        UR = col_meteorology.select("RH_G")      
        Rn24hobs = col_meteorology.select("Rn24h_G") 
        SRTM_ELEVATION = "USGS/SRTMGL1_003"
        srtm = ee.Image(SRTM_ELEVATION).clip(geometryReducer)
        z_alt = srtm.select('elevation')

        image = fexp_spec_ind(image)

        _date = ee.Date(image.get("system:time_start"))
        _hour=ee.Number(_date.get("hour"))
        _minuts=ee.Number(_date.get("minutes"))

        image = LST_DEM_correction(image, 
        z_alt,
        T_air,
        UR,
        sun_elevation,
        _hour,
        _minuts
        )
        d_cold_pixel = fexp_cold_pixel(image,geometryReducer,p_top_NDVI,p_coldest_Ts)

        n_Ts_cold = ee.Number(d_cold_pixel.get('temp'))

        image = fexp_radlong_up(image)
        image = fexp_radshort_down(image,z_alt, T_air, UR, sun_elevation)
        image = fexp_radlong_down(image, n_Ts_cold)  
        image = fexp_radbalance(image)
        image = fexp_soil_heat(image)

        d_hot_pixel=fexp_hot_pixel(image, geometryReducer, p_lowest_NDVI, p_hottest_Ts)
        date_string = _date.format("YYYY-MM-dd")

        image = fexp_sensible_heat_flux(image, ux, UR, Rn24hobs, n_Ts_cold, 
        d_hot_pixel, date_string, geometryReducer)
        image = fexp_et(image, Rn24hobs)

        # For backwards compatibility:
        image = (image.addBands(T_air)
                 .addBands(ux)
                 .addBands(UR)
                 .addBands(Rn24hobs))
        return image
    return fsebal

#IMAGE FUNCTION
class Image():

    #ENDMEMBERS DEFAULT
    #ALLEN ET AL. (2013)
    def __init__(self,
                 image,
                 NDVI_cold=5,
                 Ts_cold=20,
                 NDVI_hot=10,
                 Ts_hot=20):

        #GET INFORMATIONS FROM IMAGE
        self.image = ee.Image(image)
        self._index=self.image.get('system:index')
        self.cloud_cover=self.image.get('CLOUD_COVER')
        self.LANDSAT_ID=self.image.get('LANDSAT_ID')
        self.landsat_version= ee.String(self.image.get('SATELLITE'))
        self.azimuth_angle=self.image.get('SOLAR_ZENITH_ANGLE')
        self.time_start=self.image.get('system:time_start')
        self._date=ee.Date(self.time_start)
        self._year=ee.Number(self._date.get('year'))
        self._month=ee.Number(self._date.get('month'))
        self._day=ee.Number(self._date.get('day'))
        self._hour=ee.Number(self._date.get('hour'))
        self._minuts = ee.Number(self._date.get('minutes'))
        self.crs = self.image.projection().crs()
        self.transform = ee.List(ee.Dictionary(ee.Algorithms.Describe(self.image.projection())).get('transform'))
        self.date_string=self._date.format('YYYY-MM-dd')

        self.WRS_PATH = self.image.get("WRS_PATH")
        self.WRS_ROW = self.image.get("WRS_ROW")

        #ENDMEMBERS
        self.p_top_NDVI=ee.Number(NDVI_cold)
        self.p_coldest_Ts=ee.Number(Ts_cold)
        self.p_lowest_NDVI=ee.Number(NDVI_hot)
        self.p_hottest_Ts=ee.Number(Ts_hot)

        # GET CALIBRATED RADIANCE FROM TOA IMAGE
        # Note: no longer defining self.image_toa
        def rad_collection(img_collection):
            return (ee.ImageCollection(img_collection)
                .filterDate(self._date,self._date.advance(1,'day'))
                .filter(ee.Filter.eq("WRS_PATH", self.WRS_PATH))
                .filter(ee.Filter.eq("WRS_ROW", self.WRS_ROW))
                .filter(ee.Filter.eq("SPACECRAFT_ID", self.landsat_version))
                .map(ee.Algorithms.Landsat.calibratedRadiance)
            )
        l5_rad = ee.ImageCollection(
            rad_collection("LANDSAT/LT05/C01/T1")
        ).select([5],["T_RAD"])
        l7_rad = ee.ImageCollection(
            rad_collection("LANDSAT/LE07/C01/T1")
        ).select([5],["T_RAD"])
        l8_rad = ee.ImageCollection(
            rad_collection("LANDSAT/LC08/C01/T1")
        ).select([9],["T_RAD"])
        col_rad = l8_rad.merge(l7_rad).merge(l5_rad)
        self.col_rad = self.image.addBands(col_rad.first())

        #RENAME LANDSAT BANDS
        band_numbers = ee.Dictionary({
            "LANDSAT_5": ee.List([0,1,2,3,4,5,6,9]),
            "LANDSAT_7": ee.List([0,1,2,3,4,5,6,9]), 
            "LANDSAT_8": ee.List([0,1,2,3,4,5,6,7,10])})
        band_names = ee.Dictionary({
            "LANDSAT_5": ee.List(["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"]),
            "LANDSAT_7": ee.List(["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"]), 
            "LANDSAT_8": ee.List(["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])})

        self.image = self.image.select(
            band_numbers.get(self.landsat_version), 
            band_names.get(self.landsat_version))

        # CLOUD MASK and ALBEDO TASUMI ET AL. (2008) 
        self.image = ee.ImageCollection(
            ee.Image(
            ee.Algorithms.If(
                self.landsat_version.equals("LANDSAT_8"),
                f_albedoL8(f_cloudMaskL8_SR(self.image)),
                f_albedoL5L7(f_cloudMaskL457_SR(self.image))
            ))
        )

        sebal_algorithm = sebal(NDVI_cold=NDVI_cold, Ts_cold=Ts_cold, NDVI_hot=NDVI_hot, Ts_hot=Ts_hot)
        self.image = sebal_algorithm(self.image.first())

        # TODO -- decide whether to keep these or not in the Image object.
        # They are not really needed anymore, but keeping here for
        # backwards compatibility. However, some of these would
        # require .getInfo() to be compatible with v0.1.1.  
        self.T_air = self.image.get("AirT_G")
        self.ux = self.image.get("ux_G")
        self.UR = self.image.get("RH_G")
        self.Rn24hobs = self.image.get("Rn24h_G")
        self.geometryReducer=self.image.geometry().bounds()
        self.camada_clip=self.image.select('BRT')
        self.sun_elevation=ee.Number(90).subtract(self.azimuth_angle)
        self.srtm = ee.Image('USGS/SRTMGL1_003').clip(self.geometryReducer)
        self.z_alt = self.srtm.select('elevation')
        self.d_hot_pixel=fexp_hot_pixel(self.image, self.geometryReducer,self.p_lowest_NDVI, self.p_hottest_Ts)