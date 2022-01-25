from climada.hazard import TCTracks
import sys
import os

from config import OUT_DATA_DIR


def main(basin='EP', n_tracks=10):
    YEAR_RANGE = (1980, 2020)
    NB_SYN_TRACKS = int(n_tracks)
    # pool = Pool() # start a pathos pool
    tc_tracks = TCTracks()
    tc_tracks.read_ibtracs_netcdf(genesis_basin=basin, year_range=YEAR_RANGE)
    tc_tracks.equal_timestep(time_step_h=1)

    tc_tracks.calc_perturbed_trajectories(nb_synth_tracks=NB_SYN_TRACKS)
    tc_tracks.equal_timestep(time_step_h=1)
    path = os.path.join(OUT_DATA_DIR, "tracks", str(n_tracks), basin)
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    tc_tracks.write_netcdf(path)


if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])
