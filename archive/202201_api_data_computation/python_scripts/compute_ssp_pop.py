#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:57:40 2022

@author: kampu
"""
import os
import sys
from os import path
import numpy as np
import pandas as pd
import geopandas as gpd
from datetime import datetime
from rasterio.warp import Resampling

from climada.hazard import Centroids
from config import DATA_DIR

from climada.entity import Exposures
from climada.util.api_client import Client
from climada.util.coordinates import get_land_geometry, pts_to_raster_meta, get_resolution
from create_log_file import log_msg
from pycountry import countries
# Thanks a million to Mannie for 99% of this script

SSP_LIST = ['1', '2', '3', '4', '5']

YEAR_LIST = np.arange(2020, 2100, 20)

SSP_FILE = '/cluster/project/climate/szelie/ssp/SSP{ssp}_1km/ssp{ssp}_total_{year}.tif'

INTERMEDIATE_FILE = '/nfs/n2o/wcr/szelie/scratch/intermediate_ssp{ssp}_{year}_{country}.tif'  # intermediate files, delete after running the script

SAVE_FILE_COUNTRY = 'total_population_150as_ssp{ssp}_{year}_{country}.hdf5'

SAVE_FILE_GLOBAL = 'total_population_150as_ssp{ssp}_{year}_global.hdf5'
# RESAMPLING_METHOD = Resampling.bilinear
RESAMPLING_METHOD = Resampling.sum

CENT_FILE_PATH = os.path.join(DATA_DIR, "centroids/earth_centroids_150asland_1800asoceans_distcoast_region.hdf5")

# get all the available country list from data api
client = Client()
litpop_dataset_infos = client.list_dataset_infos(data_type='litpop')
all_properties = client.get_property_values(litpop_dataset_infos)
country_list = all_properties['country_iso3alpha']

def make_pop_ssps_countries():
    centroids_global = Centroids.from_hdf5(CENT_FILE_PATH)
    for country in country_list:
        LOG_FILE = "progress_make_ssp_pop_by_countries.txt"
        log_msg(f"Country {country} started\n", LOG_FILE)
        path = "".join([DATA_DIR, '/ssp_pop/v1/country/ssp/'])
        if np.all([
            os.path.exists(os.path.join(path, ssp, SAVE_FILE_COUNTRY.format(ssp=ssp, year=yr, country=country)))
            for ssp in SSP_LIST
            for yr in YEAR_LIST]
        ):
            log_msg(f"Country {country} already exists for all years and SSPs\n", LOG_FILE)
            continue

        # get meta data from default litpop
        try:
            # exp_litpop = client.get_litpop_default(country=country)
            country_numeric = countries.get(alpha_3=country).numeric
            centroids = centroids_global.select(reg_id=int(country_numeric))
            centroids_df = pd.DataFrame({"cent_lat": centroids.lat, "cent_lon": centroids.lon})
            centroids_df['lat_trunc'] = np.round(centroids_df['cent_lat'], 5)
            centroids_df['lon_trunc'] = np.round(centroids_df['cent_lon'], 5)

            if len(centroids.lat) > 1:
                # If all points lie on a single line of lat or lon, add another point offset so we can make a meta
                if len(np.unique(centroids.lat)) == 1 or len(np.unique(centroids.lon) == 1):
                    lat_diff = centroids.lat[1] - centroids.lat[0]
                    lon_diff = centroids.lon[1] - centroids.lon[0]
                    centroids.lat = np.append(centroids.lat, centroids.lat[-1] + lon_diff)
                    centroids.lon = np.append(centroids.lon, centroids.lon[-1] + lat_diff)

                centroids.set_lat_lon_to_meta()
                dst_meta = centroids.meta
        except Client.NoResult:
            log_msg(f"{country} not available in data api\n", LOG_FILE)
            continue

        if country == 'ATA':
            log_msg(f"{country} not available in SSP projections\n", LOG_FILE)
            continue

        # get land geometry
        country_geometry = get_land_geometry(country)

        for year in YEAR_LIST:
            log_msg(f"Year {year} begins\n", LOG_FILE)
            # if country in ['ATA', 'HMD', 'UMI', 'ESH']: # countries not available
            #     log_msg(f"{country} not available\n", LOG_FILE)
            #     continue

            if np.all([os.path.exists(
                    os.path.join(path, ssp, SAVE_FILE_COUNTRY.format(ssp=ssp, year=year, country=country))) for ssp in
                           SSP_LIST]):
                log_msg(f"Country {country} already exists for year {year} with all SSPs\n", LOG_FILE)
                continue

            for ssp in SSP_LIST:
                if os.path.exists(
                        os.path.join(path, ssp, SAVE_FILE_COUNTRY.format(ssp=ssp, year=year, country=country))):
                    log_msg(f"Country {country} SSP {ssp} already exists\n", LOG_FILE)
                    continue

                # get the population from the ssp layer reprojected
                exp_int = Exposures.from_raster(SSP_FILE.format(ssp=ssp, year=year),
                                                geometry=country_geometry)
                raw_total = np.sum(exp_int.gdf['value'])
                if len(centroids.lat) > 1:
                    # Trim centroids outside of the hazard's buffer zone. Sometimes this messes with the inferred meta
                    lat_bounds = np.sort(np.array([
                        dst_meta['transform'][5],
                        dst_meta['transform'][5] + dst_meta['transform'][4] * dst_meta['height']
                    ]))
                    lon_bounds = np.sort(np.array([
                        dst_meta['transform'][2],
                        dst_meta['transform'][2] + dst_meta['transform'][0] * dst_meta['width']
                    ]))
                    exp_int.gdf = exp_int.gdf[
                        (exp_int.gdf['latitude'] >= lat_bounds[0] - 0.00001) & \
                        (exp_int.gdf['latitude'] <= lat_bounds[1] + 0.00001) & \
                        (exp_int.gdf['longitude'] >= lon_bounds[0] - 0.00001) & \
                        (exp_int.gdf['longitude'] <= lon_bounds[1] + 0.00001)
                        ]
                    bounds = (
                        min(exp_int.gdf['longitude']),
                        min(exp_int.gdf['latitude']),
                        max(exp_int.gdf['longitude']),
                        max(exp_int.gdf['latitude'])
                    )
                    exp_int.meta['height'], exp_int.meta['width'], exp_int.meta['transform'] = pts_to_raster_meta(
                        points_bounds=bounds,
                        res=get_resolution(exp_int.gdf['longitude'], exp_int.gdf['latitude'])
                    )
                    if exp_int.meta['height'] * exp_int.meta['width'] != len(exp_int.gdf):
                        log_msg('Something is up with the meta again', LOG_FILE)
                        raise ValueError('Fix this')
                exp_int.write_raster(INTERMEDIATE_FILE.format(ssp=ssp, year=year, country=country))

                # get the exposure layer within the country boundaries
                if len(centroids.lat) == 1:
                    exp_ssp = exp_int
                    exp_ssp.meta = None
                    exp_ssp.gdf = pd.DataFrame({
                        'latitude': centroids.lat,
                        'longitude': centroids.lon,
                        'value': [np.sum(exp_ssp.gdf['value'])]
                    })
                    log_msg(f"{country} Country size has only one hazard centroid. Aggregating. \n", LOG_FILE)
                else:
                    exp_ssp = Exposures.from_raster(INTERMEDIATE_FILE.format(ssp=ssp, year=year, country=country),
                                                    dst_crs=dst_meta['crs'],
                                                    transform=dst_meta['transform'],
                                                    width=dst_meta['width'],
                                                    height=dst_meta['height'],
                                                    resampling=RESAMPLING_METHOD)
                    exp_ssp.gdf = exp_ssp.gdf[exp_ssp.gdf['value'] > 0]
                    exp_ssp.gdf['lat_trunc'] = np.round(exp_ssp.gdf['latitude'], 5)
                    exp_ssp.gdf['lon_trunc'] = np.round(exp_ssp.gdf['longitude'], 5)
                    exp_ssp.gdf = exp_ssp.gdf.merge(centroids_df, how="left", on=["lat_trunc", "lon_trunc"])
                    exp_ssp.gdf = exp_ssp.gdf[~np.isnan(exp_ssp.gdf['cent_lat'])]
                    exp_ssp.gdf = exp_ssp.gdf[~np.isnan(exp_ssp.gdf['cent_lat'])]
                    exp_ssp.gdf.drop(['lat_trunc', 'lon_trunc', 'latitude', 'longitude'], axis=1, inplace=True)
                    exp_ssp.gdf.rename({'cent_lat': 'latitude', 'cent_lon': 'longitude'}, axis=1, inplace=True)
                    regrid_total = np.sum(exp_ssp.gdf['value'])
                    regrid_change_pct = abs((regrid_total - raw_total) / raw_total) * 100
                    if regrid_change_pct > 2:
                        log_msg(
                            f"Country {country} SSP {ssp} saw a {regrid_change_pct} change in population from regridding: {raw_total} to {regrid_total}. Renormalising.\n",
                            LOG_FILE)
                    exp_ssp.gdf['value'] = exp_ssp.gdf['value'] * raw_total / regrid_total

                    exp_ssp.meta = {}
                    exp_ssp.set_crs('epsg:4326')
                    exp_ssp.gdf = gpd.GeoDataFrame(exp_ssp.gdf, geometry=gpd.points_from_xy(exp_ssp.gdf.longitude,
                                                                                            exp_ssp.gdf.latitude))
                    exp_ssp.gdf['region_id'] = countries.get(alpha_3=country).numeric
                    exp_ssp.gdf['Impf_'] = 1
                    exp_ssp.gdf = exp_ssp.gdf.dropna()
                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                exp_ssp.write_hdf5(
                    os.path.join(path, ssp, SAVE_FILE_COUNTRY.format(ssp=ssp, year=year, country=country)))
                log_msg(f"Country {country} SSP {ssp} completed\n", LOG_FILE)


def make_pop_ssps_global():
    exp_list = []
    LOG_FILE = "progress_make_ssp_pop_global.txt"
    log_msg(f"Concatenate global started\n", LOG_FILE)

    for ssp in SSP_LIST:
        for year in YEAR_LIST:
            log_msg(f"SSP{ssp} and year {year} started. \n", LOG_FILE)
            for country in country_list:
                try:
                    path = "".join([DATA_DIR, '/ssp_pop/v1/country/ssp/', ssp])
                    exp_ssp = Exposures.from_hdf5(os.path.join(path, SAVE_FILE_COUNTRY.format(ssp=ssp, year=year, country=country)))
                    exp_list.append(exp_ssp)
                except FileNotFoundError:
                    continue
            exp_global = Exposures.concat(exp_list)
            path = "".join([DATA_DIR, '/ssp_pop/v1/global/ssp/', ssp])
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            exp_global.write_hdf5(os.path.join(path,SAVE_FILE_GLOBAL.format(ssp=ssp, year=year)))


if __name__ == "__main__":
    make_pop_ssps_countries()
    make_pop_ssps_global()
