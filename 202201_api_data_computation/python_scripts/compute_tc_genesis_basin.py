from climada.hazard import Centroids, TropCyclone, TCTracks
import sys
import os
from config import DATA_DIR
from create_log_file import log_msg

basins_centroids = {'SI': (0, 140, -60, 0), 'SP': (100, -120, -60, 0),
                    'SA': (-120, 0, -60, 0), 'NA': (-120, 0, 0, 60),
                    'WP': (100, -180, 0, 60), 'NI': (0, -100, 0, 60),
                    'EP': (100, -70, 0, 60)}

CENT_FILE_PATH = os.path.join(DATA_DIR, "centroids/earth_centroids_150asland_1800asoceans_distcoast_region.hdf5")
OUT_FILE_NAME = 'tropical_cyclone_{tracks}_150arcsec_genesis_{basin}_1980_2020.hdf5'


def main(basin='EP', n_tracks=10):
    LOG_FILE = "progress_make_tc_basin.txt"
    log_msg(f"Starting computing TC for basin {basin}.\n", LOG_FILE)

    path0 = os.path.join(DATA_DIR, 'tropical_cyclones/genesis_basin', str(n_tracks) + 'synth_tracks', basin)
    centroids = Centroids()
    centroids.read_hdf5(CENT_FILE_PATH)

    tc_tracks = TCTracks()
    tc_tracks.read_netcdf(os.path.join(DATA_DIR, "tracks", str(n_tracks), basin))
    tc_haz = TropCyclone()
    centroids = centroids.select(extent=tc_tracks.get_extent(5))
    tc_haz.set_from_tracks(tc_tracks, centroids)
    tc_haz.check()
    path = os.path.join(path0, "historical")
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    tracks = str(n_tracks) + 'synth_tracks'
    file_name = OUT_FILE_NAME.format(tracks=tracks, basin=basin)
    tc_haz.write_hdf5(os.path.join(path, file_name))
    log_msg(f"Finished computing TC for basin {basin}.\n", LOG_FILE)
    tc_haz.appp

if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])