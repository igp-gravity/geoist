"""
Functions to load sample data
"""
import os
import numpy as np
import pandas as pd
import geoist
from . import datarepo
import shutil
import sys
import zipfile

Drepo = datarepo.create(
    path=["~", ".verde", "data"],
    base_url="https://raw.githubusercontent.com/igp-gravity/geoist/master/dataset/",
    version='0.0.1',
    version_dev="master",
    env="./",
	registry={"gradata.csv": "7fa07c6dd5cd56f09d9f608efd45cb2c361da38217b0725eadb1bcac39006eb1"}
)

#Drepo.load_registry(os.path.join(os.path.dirname(__file__), "registry.txt"))

Drepo1 = datarepo.create(
    path=["~", ".catalog", "data"],
    base_url="https://raw.githubusercontent.com/gravity-igpcea/dataset/master/",
    version='0.0.1',
    version_dev="master",
    env="./",
	registry={}
)

Gitee_repo_url = 'https://gitee.com/cea2020/geodataset/raw/master/demodata/'
Coding_repo_url = 'https://cea2020.coding.net/p/geodataset/d/geodataset/git/raw/master/'
ispec_catalog_url = 'http://10.2.14.222/catalog/'

def delete_downloads():
    """Delete all downloaded examples to free space or update the files."""
    shutil.rmtree(geoist.EXAMPLES_PATH)
    os.makedirs(geoist.EXAMPLES_PATH)
    return True


def _decompress(filename):
    zip_ref = zipfile.ZipFile(filename, 'r')
    zip_ref.extractall(geoist.EXAMPLES_PATH)
    return zip_ref.close()

def _get_gitee_file_url(prjname, filename):
    return 'https://gitee.com/cea2020/{}/raw/master/{}'.format(prjname, filename)

def _get_coding_file_url(prjname, dataname, filename):
    return 'https://cea2020.coding.net/p/{}/d/{}/git/raw/master/{}'.format(prjname, dataname, filename)
	
def _retrieve_file(url, filename):
    # First check if file has already been downloaded
    local_path = os.path.join(geoist.EXAMPLES_PATH, os.path.basename(filename))
    local_path_no_zip = local_path.replace('.zip', '')
    if os.path.isfile(local_path_no_zip) or os.path.isdir(local_path_no_zip):
        return local_path
    # Make sure folder exists!
    if not os.path.isdir(os.path.dirname((local_path))):
        os.makedirs(os.path.dirname((local_path)))    
	# grab the correct url retriever
    import requests
    response = requests.get(url) #, stream=True)
    # datafile = io.StringIO(response.content.decode('utf8'))
    # print(datafile)
    # print(response)
    if response.status_code == 200:
        with open(local_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size = 1024):
                f.write(chunk)
    else:
        raise IOError(response)

    if geoist.get_ext(local_path) in ['.zip']:
        _decompress(local_path)
        local_path = local_path[:-4]
    return local_path

def _download_file(filename, prjname = 'geodataset'):
    url = _get_gitee_file_url(prjname, filename)
    return _retrieve_file(url, filename)


###############################################################
def _setup_map(
    ax, xticks, yticks, crs, region, land=None, ocean=None, borders=None, states=None
):
    """
    Setup a Cartopy map with land and ocean features and proper tick labels.
    """
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

    if land is not None:
        ax.add_feature(cfeature.LAND, facecolor=land)
    if ocean is not None:
        ax.add_feature(cfeature.OCEAN, facecolor=ocean)
    if borders is not None:
        ax.add_feature(cfeature.BORDERS, linewidth=borders)
    if states is not None:
        ax.add_feature(cfeature.STATES, linewidth=states)
    ax.set_extent(region, crs=crs)
    # Set the proper ticks for a Cartopy map
    ax.set_xticks(xticks, crs=crs)
    ax.set_yticks(yticks, crs=crs)
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())


def fetch_gra_data():
    """
    Fetch sample bathymetry data from Baja California.
    """
    data_file = Drepo.fetch("gradata.csv")
    data = pd.read_csv(data_file)
    return data

def fetch_catalogCN():
    """
    Fetch Chinese earthquake catalog data from Github Repos.
    """
    
    data_file = Drepo1.fetch("cncat_mw.txt")
    data = pd.read_csv(data_file, sep='\t',encoding="utf-8")
    return data

def fetch_catalog(filename):
    """
    Fetch earthquake catalog data from Github Repos.
    """
    
    data_file = Drepo1.fetch(filename)
    data = pd.read_csv(data_file, sep='\t',encoding="utf-8")
    return data

def fetch_catalogGEM():
    """
    Fetch GEM Global earthquake catalog data from Github Repos.
    """
    data_file = Drepo1.fetch("isc-gem-cat.txt")
    data = pd.read_csv(data_file, sep='\t',encoding="utf-8")
    return data


def setup_cn_bathymetry_map(
    ax, region=(100.0, 112.0, 40.0, 55.0), land="gray", ocean=None
):
    """
    Setup a Cartopy map for the bathymetry dataset.

    Parameters
    ----------
    ax : matplotlib Axes
        The axes where the map is being plotted.
    region : list = [W, E, S, N]
        The boundaries of the map region in the coordinate system of the data.
    land : str or None
        The name of the color of the land feature or None to omit it.
    ocean : str or None
        The name of the color of the ocean feature or None to omit it.

    """
    import cartopy.crs as ccrs

    _setup_map(
        ax,
        xticks=np.arange(100, 111, 5),
        yticks=np.arange(45, 55, 2),
        land=land,
        ocean=ocean,
        region=region,
        crs=ccrs.PlateCarree(),
    )



def fetch_baja_bathymetry():
    data_file = Drepo.fetch("baja-bathymetry.csv.xz")
    data = pd.read_csv(data_file, compression="xz")
    return data


def setup_baja_bathymetry_map(
    ax, region=(245.0, 254.705, 20.0, 29.99), land="gray", ocean=None
):
    import cartopy.crs as ccrs

    _setup_map(
        ax,
        xticks=np.arange(-114, -105, 2),
        yticks=np.arange(21, 30, 2),
        land=land,
        ocean=ocean,
        region=region,
        crs=ccrs.PlateCarree(),
    )


def fetch_rio_magnetic():
    data_file = Drepo.fetch("rio-magnetic.csv.xz")
    data = pd.read_csv(data_file, compression="xz")
    return data


def setup_rio_magnetic_map(ax, region=(-42.6, -42, -22.5, -22)):
    import cartopy.crs as ccrs

    _setup_map(
        ax,
        xticks=np.arange(-42.5, -42, 0.1),
        yticks=np.arange(-22.5, -21.99, 0.1),
        land=None,
        ocean=None,
        region=region,
        crs=ccrs.PlateCarree(),
    )


def fetch_california_gps():
    data_file = Drepo.fetch("california-gps.csv.xz")
    data = pd.read_csv(data_file, compression="xz")
    return data


def setup_california_gps_map(
    ax, region=(235.2, 245.3, 31.9, 42.3), land="gray", ocean="skyblue"
):
    import cartopy.crs as ccrs

    _setup_map(
        ax,
        xticks=np.arange(-124, -115, 4),
        yticks=np.arange(33, 42, 2),
        land=land,
        ocean=ocean,
        region=region,
        crs=ccrs.PlateCarree(),
    )


def fetch_texas_wind():
    data_file = Drepo.fetch("texas-wind.csv")
    data = pd.read_csv(data_file)
    return data


def setup_texas_wind_map(
    ax, region=(-107, -93, 25.5, 37), land="#dddddd", borders=0.5, states=0.1
):

    import cartopy.crs as ccrs

    _setup_map(
        ax,
        xticks=np.arange(-106, -92, 3),
        yticks=np.arange(27, 38, 3),
        land=land,
        ocean=None,
        region=region,
        borders=borders,
        states=states,
        crs=ccrs.PlateCarree(),
    )
