from climada.hazard import TCTracks
import sys
import os
from config import DATA_DIR


def main(basin='EP', n_tracks=10):
    year_range = (1980, 2020)
    nb_syn_tracks = int(n_tracks)
    tc_tracks = TCTracks()
    tc_tracks.read_ibtracs_netcdf(genesis_basin=basin, year_range=year_range)
    tc_tracks.equal_timestep(time_step_h=1)
    if nb_syn_tracks>0:
        tc_tracks.calc_perturbed_trajectories(nb_synth_tracks=nb_syn_tracks)
    path = os.path.join(DATA_DIR, "tracks", str(n_tracks), basin)
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    tc_tracks.write_netcdf(path)


if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])
