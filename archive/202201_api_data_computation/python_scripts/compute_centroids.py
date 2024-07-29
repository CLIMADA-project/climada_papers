import os
from climada.hazard import Centroids
import shapely
from shapely.ops import unary_union
import shapely.geometry
import cartopy.io.shapereader as shpreader
from create_log_file import log_msg


def make_base_centroids(out_file_path, bounds=(-180, -60, 180, 60), res_land_arcsec=150, res_ocean_arcsec=1800,
                        land_buffer=0.1):
    """
    Create centroids for a given bounding box, with given resolution for land and ocean, and then combine both.

    The function reads the geographical boundaries, creates centroids for both land and ocean,
    combines them, computes additional attributes, and saves the centroids to an output file.

    Parameters:
        out_file_path (str): Path for saving the output centroids.
        bounds (tuple, optional): Boundaries for creating centroids in the format (lon_min, lat_min, lon_max, lat_max).
                                  Default is (-180, -60, 180, 60).
        res_land_arcsec (int, optional): Resolution for land in arcseconds. Default is 150.
        res_ocean_arcsec (int, optional): Resolution for ocean in arcseconds. Default is 1800.
        land_buffer (float, optional): Buffer around land in degrees, to still consider the land resolution.
        Default is 0.1.

    Returns:
        None
    """
    LOG_FILE = "progress_make_centroids.txt"
    log_msg(f"Started with bounds {bounds}, resolution on land {res_land_arcsec},"
            f"and resolution on ocean {res_ocean_arcsec}\n", LOG_FILE)

    res_land = res_land_arcsec / 3600 #getting resolution in degrees
    res_ocean = res_ocean_arcsec / 3600
    cent_land = Centroids.from_pnt_bounds(points_bounds=bounds, res=res_land)
    cent_land.set_meta_to_lat_lon()
    cent_ocean = Centroids.from_pnt_bounds(points_bounds=bounds, res=res_ocean)
    cent_ocean.set_meta_to_lat_lon()
    shpfilename = shpreader.natural_earth(category='physical', name='land', resolution='10m')
    land = shpreader.Reader(shpfilename)
    land_geometrys = [x.geometry for x in land.records()]
    cu = unary_union(list(land_geometrys))
    cu = cu.buffer(land_buffer, resolution=10) # creating a 10km buffer around land

    #making an ocean layer
    mask = shapely.vectorized.contains(cu, cent_ocean.lon, cent_ocean.lat)
    cent_ocean.lat = cent_ocean.lat[~mask]
    cent_ocean.lon = cent_ocean.lon[~mask]

    #making a land layer
    mask = shapely.vectorized.contains(cu, cent_land.lon, cent_land.lat)
    cent_land.lat = cent_land.lat[mask]
    cent_land.lon = cent_land.lon[mask]

    cent = cent_land
    cent.append(cent_ocean) #combining both
    log_msg(f"centroids ready, starting computation of other attributes", LOG_FILE)
    cent.set_region_id()
    log_msg(f"region if set, starting set_on_land()", LOG_FILE)

    cent.set_on_land()
    cent = cent.select(extent=(bounds[0], bounds[2], bounds[1], bounds[3]))
    log_msg(f"set_on_land done, setting dist to coast", LOG_FILE)

    cent.set_dist_coast() #takes a lot of time
    cent.check()
    cent.write_hdf5(out_file_path)
    log_msg(f"done, file saved", LOG_FILE)


file_name = os.path.join('centroids', 'earth_centroids_150asland_1800asoceans_distcoast_region.hdf5')
file_name_no_poles = os.path.join('centroids',
                          'earth_centroids_150asland_1800asoceans_distcoast_region_nopoles.hdf5')

if __name__ == "__main__":
    make_base_centroids(file_name, bounds=(-180, -90, 180, 90))
    make_base_centroids(file_name_no_poles, bounds=(-180, -60, 180, 60))
