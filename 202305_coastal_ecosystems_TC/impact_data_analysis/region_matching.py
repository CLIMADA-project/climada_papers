"""
Script to match impact data with shapefile containing regions
Author: Sarah Hülsen
"""

from pathlib import Path
import pandas as pd
import geopandas as gpd

years = ['2000', '2020']
proj_eck4 = '+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'  # World Eckert IV


for year in years:
    path = Path('../results/final/')
    region_path = Path('../data/')
    # load geodataframe containing regions
    reg = gpd.read_file(f'{region_path}UN_regions_custom.gpkg')

    # load csv containing all exposure points
    imp = pd.read_csv(f'{path}all_basins_hist_STORM_{year}_wp.csv')

    # convert exposure df to gdf
    imp_gdf = gpd.GeoDataFrame(imp, geometry=gpd.points_from_xy(imp.longitude, imp.latitude))
    imp_gdf = imp_gdf.set_crs("EPSG:4326")

    # project geodataframes
    imp_gdf = imp_gdf.to_crs(proj_eck4)
    reg = reg.to_crs(proj_eck4)

    # perform spatial join
    sjoin = gpd.sjoin_nearest(left_df=imp_gdf, right_df=reg, how='left', max_distance=100000, distance_col='distance_m')
    sjoin = sjoin.to_crs("EPSG:4326")

    # create analysis file for merging with other impact files
    # select only necessary columns
    cols = ['longitude',
            'latitude',
            'ISO3',
            'UN',
            'NAME',
            'REGION',
            'SUBREGION',
            'sub_region_name',
            'region_name',
            'custom_region',
            'distance_m']
    df = sjoin[cols]

    # rename columns
    df.columns = ['longitude',
                      'latitude',
                      'ISO3',
                      'UN_country_code',
                      'country_name',
                      'UN_region_code',
                      'UN_subregion_code',
                      'sub_region_name',
                      'region_name',
                      'custom_region',
                      'distance_m']

    # save file
    df.to_csv(f'{path}custom_regions_analysis_file_{year}.csv')
