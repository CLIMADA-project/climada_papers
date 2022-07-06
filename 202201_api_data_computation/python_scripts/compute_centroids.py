import os
from climada.hazard import Centroids
import numpy as np
import geopandas as gpd
import shapely
from shapely.ops import unary_union
import shapely.geometry
import cartopy.io.shapereader as shpreader
from config import OUT_DATA_DIR


def make_base_centroids(out_file_path, bounds=(-180, -60, 180, 60), res_land_arcsec=150, res_ocean_arcsec=1800,
                        land_buffer=0.1):
    res_land = res_land_arcsec/3600
    res_ocean = res_ocean_arcsec/3600
    cent_land = Centroids.from_pnt_bounds(points_bounds=bounds, res=res_land)
    cent_land.set_meta_to_lat_lon()
    cent_ocean = Centroids.from_pnt_bounds(points_bounds=bounds, res=res_ocean)
    cent_ocean.set_meta_to_lat_lon()
    shpfilename = shpreader.natural_earth(category='physical', name='land', resolution='10m')
    land = shpreader.Reader(shpfilename)
    land_geometrys = [x.geometry for x in land.records()]
    cu = unary_union(list(land_geometrys))
    cu = cu.buffer(land_buffer, resolution=10)

    mask = shapely.vectorized.contains(cu, cent_ocean.lon, cent_ocean.lat)
    cent_ocean.lat = cent_ocean.lat[~mask]
    cent_ocean.lon = cent_ocean.lon[~mask]

    mask = shapely.vectorized.contains(cu, cent_land.lon, cent_land.lat)
    cent_land.lat = cent_land.lat[mask]
    cent_land.lon = cent_land.lon[mask]

    cent = cent_land
    cent.append(cent_ocean)
    cent.set_region_id()
    #cent.set_geometry_points()
    #cent.set_lat_lon_to_meta()
    cent.set_on_land()
    cent = cent.select(extent=(bounds[0], bounds[2], bounds[1], bounds[3]))
    cent.set_dist_coast(scheduler="threads")
    #cent.set_geometry_points()
    cent.check()
    cent.write_hdf5(out_file_path)



file_name = os.path.join(OUT_DATA_DIR, 'centroids','earth_centroids_150asland_1800asoceans_distcoast_region2.hdf5')
if __name__ == "__main__":
    make_base_centroids(file_name, bounds=(-180, -90, 180, 90))


file_name = os.path.join(OUT_DATA_DIR, 'centroids','earth_centroids_150asland_1800asoceans_distcoast_region_nopoles2.hdf5')
if __name__ == "__main__":
    make_base_centroids(file_name, bounds=(-180, -60, 180, 60))