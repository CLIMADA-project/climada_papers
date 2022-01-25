from climada.hazard import Centroids, TropCyclone, TCTracks
import sys
import os
from config import OUT_DATA_DIR


def main(basin='EP', climate_scenarios=[26, 60, 45, 85], future_years=[2040, 2060, 2080], n_tracks=10):
    path0 = os.path.join(OUT_DATA_DIR, '/tropical_cyclones/genesis_basin', str(n_tracks) + 'synth_tracks', basin)
    centroids = Centroids()
    centroids.read_hdf5(os.path.join(OUT_DATA_DIR,
                                     "centroids/earth_centroids_150asland_1800water_distcoast_region_nopoles.hdf5"))

    #if basin in ("SI", "SP", "SA"):
    #    centroids = centroids.select(extent=(-180, 180, -71, 0))
    #else:
    #    centroids = centroids.select(extent=(-180, 180, 0, 61))

    #pool = Pool()  # start a pathos pool

    tc_tracks = TCTracks()
    tc_tracks.read_netcdf(os.path.join("../data/tracks/", str(n_tracks), basin))
    tc_haz = TropCyclone()
    #centroids_cut = centroids.select(extent=tc_tracks.get_extent(5))
    tc_haz.set_from_tracks(tc_tracks, centroids)
    tc_haz.check()
    path = os.path.join(path0, "historical")
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    tracks = str(n_tracks) + 'synth_tracks'
    file_name = "_".join(['tropical_cyclone', tracks, '150arcsec_genesis', basin, '1980_2020.hdf5'])
    tc_haz = TropCyclone.from_hdf5(os.path.join(path, file_name))
    for climate_scenario in climate_scenarios:
        for year in future_years:
            rcp_str = 'rcp' + str(climate_scenario)
            path = os.path.join(path0, rcp_str, str(year))
            file_name = "_".join(['tropical_cyclone', tracks, '150arcsec', rcp_str, 'genesis', basin, str(year)])
            file_name = ".".join([file_name, 'hdf5'])
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            tc_haz_future = tc_haz.set_climate_scenario_knu(ref_year=year, rcp_scenario=climate_scenario)
            tc_haz_future.write_hdf5(os.path.join(path, file_name))
    # pool.close()
    # pool.join()


if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])