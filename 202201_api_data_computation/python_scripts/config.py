import platform

if platform.system() == 'Darwin':
    OUT_DATA_DIR = "/Users/szelie/OneDrive - ETH Zurich/data/climada_api"
    IN_DATA_DIR = "/Users/szelie/OneDrive - ETH Zurich/data/climada_api"

else:
    OUT_DATA_DIR = "/cluster/project/climate/szelie/CLIMADA_api_data"
    IN_DATA_DIR = "/cluster/project/climate/szelie/CLIMADA_api_data"