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

#SURFACE REFLECTANCE
#ATMOSPHERICALLY CORRECTED

def set_landsat_index(img):
      """Keeps system:index as LANDSAT_INDEX so that it doesn't get lost
       by merging or joining collections"""
      return img.set({"LANDSAT_INDEX": img.get("system:index")})

# Filter collection by path, row, coordinate, and cloud cover
def fexp_collection_filter(collection, start_date, end_date, th_cloud_cover, 
    n_path=None, n_row=None, coordinate = None):
    collection = (collection.filterDate(start_date, end_date)
            .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover)
    )
    if n_path: collection = collection.filterMetadata('WRS_PATH', 'equals', n_path)
    if n_row: collection = collection.filterMetadata('WRS_ROW', 'equals', n_row)
    if coordinate: collection = collection.filterBounds(coordinate)
    return collection

# GET T_RAD band from T1 products  
def fexp_trad_8(start_date,end_date,th_cloud_cover, **kwargs):
    return (fexp_collection_filter(
        ee.ImageCollection('LANDSAT/LC08/C01/T1'),
        start_date, end_date, th_cloud_cover, **kwargs)
            .map(ee.Algorithms.Landsat.calibratedRadiance)
            .map(set_landsat_index)
            .select([9],["T_RAD"]))

def fexp_trad_7(start_date,end_date,th_cloud_cover, **kwargs):
    return (fexp_collection_filter(
        ee.ImageCollection('LANDSAT/LE07/C01/T1'),
        start_date, end_date, th_cloud_cover, **kwargs)
            .map(ee.Algorithms.Landsat.calibratedRadiance)
            .map(set_landsat_index)
            .select([5],["T_RAD"]))

def fexp_trad_5(start_date,end_date,th_cloud_cover, **kwargs):
    return (fexp_collection_filter(
        ee.ImageCollection('LANDSAT/LT05/C01/T1'),
        start_date, end_date, th_cloud_cover, **kwargs)
            .map(ee.Algorithms.Landsat.calibratedRadiance)
            .map(set_landsat_index)
            .select([5],["T_RAD"]))

#GET LANDSAT 8 COLLECTION
def fexp_landsat_8(start_date,end_date, th_cloud_cover, **kwargs):
    return (fexp_collection_filter(
        ee.ImageCollection('LANDSAT/LC08/C01/T1_SR'),
        start_date, end_date, th_cloud_cover, **kwargs)
        .select([0,1,2,3,4,5,6,7,10],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
        .map(set_landsat_index))

#GET LANDSAT 7 COLLECTION
def fexp_landsat_7(start_date,end_date, th_cloud_cover, **kwargs):
    return (fexp_collection_filter(
        ee.ImageCollection('LANDSAT/LE07/C01/T1_SR'),
        start_date, end_date, th_cloud_cover, **kwargs)
        .select([0,1,2,3,4,5,6,9], ["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"])
        .map(set_landsat_index))

#GET LANDSAT 5 COLLECTION
def fexp_landsat_5(start_date,end_date, th_cloud_cover, **kwargs):
    return (fexp_collection_filter(
        ee.ImageCollection('LANDSAT/LT05/C01/T1_SR'),
        start_date, end_date, th_cloud_cover, **kwargs)
        .select([0,1,2,3,4,5,6,9], ["B","GR","R","NIR","SWIR_1","BRT","SWIR_2", "pixel_qa"])
        .map(set_landsat_index))
