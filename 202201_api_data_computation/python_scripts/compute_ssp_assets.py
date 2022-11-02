#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 17:52:10 2022

@author: evelynm
"""
import copy

import pandas as pd
import logging
from config import DATA_DIR
import os
import climada.util.coordinates as u_coord
from climada.entity.exposures import Exposures
from climada.util.constants import SYSTEM_DIR
from pycountry import countries
import numpy as np
from create_log_file import log_msg

LOGGER = logging.getLogger(__name__)
LOG_FILE = "progress_make_litpop_ssp.txt"
FILE_LITPOP = 'LitPop_150arcsec_{country}.hdf5'
SAVE_FILE_COUNTRY = 'LitPop_150as_ssp{ssp}_{year}_{country}.hdf5'
SAVE_FILE_GLOBAL = 'LitPop_150as_ssp{ssp}_{year}_global.hdf5'
FILE_SSP_MULT = os.path.join(DATA_DIR, 'ssp_assets',
                             'ssps_gdp_annual.csv')

DIR = os.path.join(DATA_DIR, 'litpop', 'countries')
SSP_LIST = ['1', '2', '3', '4', '5']

YEAR_LIST = np.arange(2020, 2100, 20)


def get_gdp_scen(country, year, ssp, file_path_factor, model='IIASA'):
    """
    Lookup function for a country's growth factor in year X, relative to the
    base year 2020, according to an SSP scenario and a modelling source.

    Annual growth factors were calculated from the SSP public database (v2.0)
    Keywan Riahi, Detlef P. van Vuuren, Elmar Kriegler, Jae Edmonds,
    Brian C. O’Neill, Shinichiro Fujimori, Nico Bauer, Katherine Calvin,
    Rob Dellink, Oliver Fricko, Wolfgang Lutz, Alexander Popp,
    Jesus Crespo Cuaresma, Samir KC, Marian Leimbach, Leiwen Jiang, Tom Kram,
    Shilpa Rao, Johannes Emmerling, Kristie Ebi, Tomoko Hasegawa, Petr Havlík,
    Florian Humpenöder, Lara Aleluia Da Silva, Steve Smith, Elke Stehfest,
    Valentina Bosetti, Jiyong Eom, David Gernaat, Toshihiko Masui, Joeri Rogelj,
    Jessica Strefler, Laurent Drouet, Volker Krey, Gunnar Luderer, Mathijs Harmsen,
    Kiyoshi Takahashi, Lavinia Baumstark, Jonathan C. Doelman, Mikiko Kainuma,
    Zbigniew Klimont, Giacomo Marangoni, Hermann Lotze-Campen, Michael Obersteiner,
    Andrzej Tabeau, Massimo Tavoni.
    The Shared Socioeconomic Pathways and their energy, land use, and
    greenhouse gas emissions implications: An overview, Global Environmental
    Change, Volume 42, Pages 153-168, 2017,
    ISSN 0959-3780, DOI:110.1016/j.gloenvcha.2016.05.009
    Selection: 1. Region - all countries, 2. Scenarios - GDP (PIK, IIASA, OECD),
    3. Variable - GDP (growth PPP)

    Parameters
    ----------
    country : str
        iso3alpha (e.g. 'JPN'), or English name (e.g. 'Switzerland')
    year : int
        The yer for which to get a GDP projection for. Any among [2020, 2099].
    ssp : int
        The SSP scenario for which to get a GDP projecion for. Any amon [1, 5].
    model : str
        The model source from which the GDP projections have been calculated.
        Either IIASA, PIK or OECD. Default is IIASA.

    Returns
    -------
    float
        The country's GDP growth relative to the year 2020, according to chosen
        SSP scenario and source.

    Example
    -------
    get_gdp_scen('Switzerland', 2067, 2)
    """

    model_long = 'IIASA GDP'
    if model == 'OECD':
        model_long = 'OECD Env-Growth'
    elif model == 'PIK':
        model_long = 'PIK GDP-32'

    ssp_long = f'SSP{str(ssp)}'

    try:
        iso3a = u_coord.country_to_iso(country, representation="alpha3")
    except LookupError:
        LOGGER.error('Country not identified: %s.', country)
        return None

    df_csv = pd.read_csv(file_path_factor,
                         usecols=['Model', 'Scenario', 'Region', str(year)])

    sel_bool = ((df_csv.Model == model_long) &
                (df_csv.Scenario == ssp_long) &
                (df_csv.Region == iso3a))
    if np.sum(sel_bool) == 0:
        return

    return df_csv.loc[sel_bool][str(year)].values[0]


def make_ssp_assets():
    LOG_FILE = "progress_make_ssp_assets.txt"
    exposures_list = []
    missing_coutry = []
    for year in YEAR_LIST:
        for ssp in SSP_LIST:
            path = os.path.join(DATA_DIR, 'ssp_assets', 'global', "".join(['ssp', ssp]))
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            file_path_global = os.path.join(path, SAVE_FILE_GLOBAL.format(ssp=ssp, year=year))
            isExist = os.path.exists(file_path_global)
            if isExist:
                log_msg(f"All countries already computed for SSP{ssp} and year {year}, continuing. \n", LOG_FILE)
                continue
            for country in countries:
                country = country.alpha_3
                log_msg(f"Starting country {country} for SSP{ssp} and year {year}. \n", LOG_FILE)
                path = os.path.join(DATA_DIR, 'ssp_assets', 'country', "".join(['ssp', ssp]))

                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                file_path = os.path.join(path, SAVE_FILE_COUNTRY.format(ssp=ssp, year=year, country=country))
                file = FILE_LITPOP.format(country=country)
                path_litpop = os.path.join(DIR, 'default', file)
                try:
                    litpop_present = Exposures.from_hdf5(path_litpop)
                except FileNotFoundError:
                    missing_coutry.append(country)
                    log_msg(f"LitPop for country {country} not fount. \n", LOG_FILE)
                    continue
                file_path_factor = FILE_SSP_MULT
                factor = get_gdp_scen(country, year, ssp, file_path_factor)
                if factor is None:
                    missing_coutry.append(country)
                    log_msg(f"SSP scenarios for {country} were not found. \n", LOG_FILE)
                    continue
                litpop_future = copy.deepcopy(litpop_present)
                litpop_future.gdf.value = litpop_future.gdf.value*factor
                log_msg(f"Data directory is  {DATA_DIR}\n", LOG_FILE)


                log_msg(f"Saving file as {file_path}\n", LOG_FILE)

                litpop_future.write_hdf5(file_path)
                log_msg(f"Finished country {country}\n", LOG_FILE)
                exposures_list.append(litpop_future)
            log_msg(f"Done creating country files, making global file. \n", LOG_FILE)
            log_msg(f"The following countries could not be created {missing_coutry}. \n", LOG_FILE)

            litpop_future_global = Exposures.concat(exposures_list)
            litpop_future_global.write_hdf5(file_path_global)
            log_msg(f"Done. \n", LOG_FILE)


if __name__ == "__main__":
    make_ssp_assets()