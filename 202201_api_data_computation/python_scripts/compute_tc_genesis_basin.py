from climada.hazard import Centroids, TropCyclone, TCTracks
import sys
import os
from config import DATA_DIR
from create_log_file import log_msg
from climada.util.api_client import Client

basins_centroids = {'SI': (0, 140, -60, 0), 'SP': (100, -120, -60, 0),
                    'SA': (-120, 0, -60, 0), 'NA': (-120, 0, 0, 60),
                    'WP': (100, -180, 0, 60), 'NI': (0, -100, 0, 60),
                    'EP': (100, -70, 0, 60)}

OUT_FILE_NAME = 'tropical_cyclone_{tracks}_150arcsec_genesis_{basin}_1980_2020.hdf5'

def main(basin='EP', n_tracks=10):
    """
    Create a TropCyclone hazard set based on TCTracks for a specific basin. Should be ran
    after compute_tc_tracks.py

    The function reads the TC tracks for the specified basin, calculates the hazard,
    and saves the TropCyclone hazard set to an output file.

    Parameters:
        basin (str, optional): Name of the TC basin. Default is 'EP'.
        n_tracks (int, optional): Number of synthetic tracks to consider. Default is 10.

    """

    client = Client()
    LOG_FILE = "progress_make_tc_basin.txt"
    log_msg(f"Starting computing TC for basin {basin}.\n", LOG_FILE)

    path0 = os.path.join(DATA_DIR, 'tropical_cyclones_v3/genesis_basin', str(n_tracks) + 'synth_tracks', basin)
    centroids = client.get_centroids(res_arcsec_land=150,
    res_arcsec_ocean=1800,
    extent=(-180, 180, -90, 90),)

    tc_tracks = TCTracks.from_netcdf(os.path.join(DATA_DIR, "tracks", str(n_tracks), basin))
    centroids = centroids.select(extent=tc_tracks.get_extent(5))
    tc_haz = TropCyclone.from_tracks(tc_tracks, centroids, pool=None)
    tc_haz.check()
    path = os.path.join(path0, "historical")
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    tracks = str(n_tracks) + 'synth_tracks'
    file_name = OUT_FILE_NAME.format(tracks=tracks, basin=basin)
    tc_haz.write_hdf5(os.path.join(path, file_name))
    log_msg(f"Finished computing TC for basin {basin}.\n", LOG_FILE)

if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])