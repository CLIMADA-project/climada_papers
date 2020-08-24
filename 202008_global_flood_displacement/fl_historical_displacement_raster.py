#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the CLIMADA papers repository:
    https://github.com/CLIMADA-project/climada_papers

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Created on Mon Aug 24 09:56:35 2020

Description:
    This script is to compute the ABSOLUTE number of displaced people due to river floods at 5km resolution globally.
    The output raster files represent the number of displaced people at every 10 years from 1970-2000 using the baseline population at year 2000 and historical simulation of flooding hazards.
    The flow of the script is similar to fl_cc_displacement_raster.py
    Please refer to the paper (submitted) and README for more information.
    For the baseline calculation please refer to flood_historical_displacement_raster.py.
    JUPYTER NOTEBOOK documented the calculation using Vietnam as a case study.

@author: manniekam
"""
