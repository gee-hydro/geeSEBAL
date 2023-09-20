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
from datetime import date

#FOLDERS
from .landsatcollection import (fexp_landsat_5,fexp_landsat_7, fexp_landsat_8,
fexp_trad_5, fexp_trad_7, fexp_trad_8)
from .masks import (f_cloudMaskL457_SR,f_cloudMaskL8_SR,f_albedoL5L7,f_albedoL8)
from .image import sebal

#COLLECTION FUNCTION
class Collection():

    #ENDMEMBERS DEFAULT
    #ALLEN ET AL. (2013)
    def __init__(self,
                 year_i,
                 month_i,
                 day_i,
                 year_e,
                 month_e,
                 day_e,
                 cloud_cover,
                 path=None,
                 row=None,
                 coordinate=None,
                 NDVI_cold=5,
                 Ts_cold=20,
                 NDVI_hot=10,
                 Ts_hot=20, 
                 max_iterations=15):

        if not any([all([path,row]),coordinate]):
            raise ValueError("Either coordinate or a combination of path/row should be specified.")
        # ^ Otherwise user could request too many Landsat images by accident. 
        spatial_filters = {"n_path":path, "n_row":row, "coordinate":coordinate}

        #INFORMATIONS
        self.path=path
        self.row=row
        self.coordiate=coordinate
        self.cloud_cover=cloud_cover
        self.start_date = ee.Date.fromYMD(year_i,month_i,day_i)
        self.i_date=date(year_i,month_i,day_i)
        self.end_date=date(year_e,month_e,day_e)
        self.n_search_days=self.end_date - self.i_date
        self.n_search_days=self.n_search_days.days
        self.end_date = self.start_date.advance(self.n_search_days, 'day')
        self.max_iterations = max_iterations

        #COLLECTIONS 
        self.collection_l5=fexp_landsat_5(self.start_date, self.end_date, self.cloud_cover, **spatial_filters)
        self.collection_l7=fexp_landsat_7(self.start_date, self.end_date, self.cloud_cover, **spatial_filters)
        self.collection_l8=fexp_landsat_8(self.start_date, self.end_date, self.cloud_cover, **spatial_filters)
        rad_l5 = fexp_trad_5(self.start_date, self.end_date, self.cloud_cover, **spatial_filters)
        rad_l7 = fexp_trad_7(self.start_date, self.end_date, self.cloud_cover, **spatial_filters)
        rad_l8 = fexp_trad_8(self.start_date, self.end_date, self.cloud_cover, **spatial_filters)
        rad_collection = rad_l5.merge(rad_l7).merge(rad_l8)

        #LIST OF IMAGES
        self.sceneListL5 = self.collection_l5.aggregate_array('system:index')
        self.sceneListL7 = self.collection_l7.aggregate_array('system:index')
        self.sceneListL8 = self.collection_l8.aggregate_array('system:index')

        #ALBEDO AND MASKS    
        self.collection_l5 = self.collection_l5.map(f_albedoL5L7).map(f_cloudMaskL457_SR)   
        self.collection_l7 = self.collection_l7.map(f_albedoL5L7).map(f_cloudMaskL457_SR)
        self.collection_l8= self.collection_l8.map(f_albedoL8).map(f_cloudMaskL8_SR)

        #JOIN COLLECTIONS
        self.collection = self.collection_l5.merge(self.collection_l7).merge(self.collection_l8)  
        filter_landsat_index =  ee.Filter.equals(leftField='LANDSAT_INDEX', rightField='LANDSAT_INDEX') 
        join_by_landsat_index = ee.ImageCollection(
            ee.Join.inner().apply(
                self.collection, 
                rad_collection, 
                filter_landsat_index
            ))
        def get_img(feature):
            return ee.Image.cat(feature.get("primary"), feature.get("secondary"))
        
        self.collection = join_by_landsat_index.map(get_img)  
        self.count = self.collection.size()
        self.CollectionList=self.collection.aggregate_array('system:index')
        self.CollectionList_image = self.collection.aggregate_array('system:index')

        # Map the sebal algorithm to the image collection
        sebal_algorithm = sebal(NDVI_cold=NDVI_cold, Ts_cold=Ts_cold, 
                                NDVI_hot=NDVI_hot, Ts_hot = Ts_hot,
                                max_iterations=max_iterations)
        self.collection = self.collection.map(sebal_algorithm).sort("system:time_start")

        # Backwards compatibility: Collection_ET is an ee.Image
        # with only the ET_24h band but renamed to {LANDSAT_INDEX}:
        self.Collection_ET = (self.collection.select("ET_24h").toBands()
                              .rename(
                                  self.collection.aggregate_array("LANDSAT_INDEX")
                              ))