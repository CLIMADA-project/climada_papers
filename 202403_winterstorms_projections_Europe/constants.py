#### define constants to be use in the analysis

##imports
#from func_preproc_climada import *
### naming
##variables:
#day:surface wind max = SWM, surface wind (mean) = SW
#Amon: near surface air temperature: TAS, zonal wind: UA, air temperature: TA
#Omon: sea-surface temperature: TOS
#processing: 98th quantile = q98, gust factor 1.67 = gst1-67, NA20 = "threshold=20" in dropna(), cal1 = type 1 calibration
#space resolution: original model resolution = br, interpolated 25deg = itp2_5
#time res
#domain: Europe = EU
#season: extended winter = winE
#base name save file = metvar + processings + time res + space res + domain + season
## climada
#CLIMADA filetypes: hazard = haz, exposures = exp, impacts = imp
#impact function: see shortnames

### constants
#aai_agg_emdat = 1000*2165734.375
aai_agg_emdat = 1000*3073582.5
CubEOT_corr_fact = 0.0015956710965959711

## model sets
modlist_allscen = ['CanESM5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3-Veg','EC-Earth3-Veg-LR','IPSL-CM6A-LR','MIROC-ES2L',
           'UKESM1-0-LL','MRI-ESM2-0','FGOALS-g3','ACCESS-ESM1-5','MIROC6','MPI-ESM1-2-LR','KACE-1-0-G'] #complete model list

modlist_ssp585 = ['AWI-CM-1-1-MR','BCC-CSM2-MR','CNRM-CM6-1-HR','EC-Earth3','EC-Earth3-CC','HadGEM3-GC31-LL','GISS-E2-1-G','GFDL-CM4','CMCC-CM2-SR5','CMCC-ESM2','HadGEM3-GC31-MM','NESM3','MPI-ESM1-2-HR','INM-CM4-8','INM-CM5-0','ACCESS-CM2']

modlist_1cen = ['CanESM5','CNRM-CM6-1','EC-Earth3-Veg','IPSL-CM6A-LR','UKESM1-0-LL','MRI-ESM2-0','FGOALS-g3','ACCESS-ESM1-5','MIROC6','MPI-ESM1-2-LR','AWI-CM-1-1-MR','BCC-CSM2-MR','CNRM-CM6-1-HR','GISS-E2-1-G','GFDL-CM4','CMCC-ESM2','HadGEM3-GC31-MM','NESM3','INM-CM5-0']

#scenarios
scenlist= ['historical','ssp126','ssp245','ssp370','ssp585']

#region names
reglist3 = ['BI','IP','FR','WEU','MED','SC','EEU']


#impact function
impflist = ["CubEOT","Em2011","Wk2021"]

## paths
#in
pathcmip6 = '/home/lseverino/MT/cmip6/' #path of the cmip6 data
pathera5 = '/home/lseverino/MT/era5/'
pathcirc = '/home/lseverino/MT/circulation/'
pathin = '/home/lseverino/MT/inputdata/'
#out
#pathcal = '/home/lseverino/MT/scripts/calibration/'
#pathfig = '/home/lseverino/MT/scripts/results/figures'
#pathhaz = "/home/lseverino/MT/scripts/results/hazards/"
#pathimp = "/home/lseverino/MT/scripts/results/impacts/"

pathcal = './calibration/'
pathfig = './results/figures'
pathhaz = "./results/hazards/"
pathimp = "./results/impacts/"

#dict for abbreviations of the cmip6 variables names
cmip6vars = {'sfcWindmax':'SWM','sfcWind':'SW','psl':'SLP','tas':'TAS','ua':'UA','ta':'TA','tos':'TOS'}

#dic for impf shortnames
impf_sht_names = {'Cubic excess-over-threshold':'CubEOT','Scaled sigmoid':'ScSig','Emanuel 2011':'Em2011','Welker 2021':'Wk2021','Schwierz 2010' : 'Sw2010'}



