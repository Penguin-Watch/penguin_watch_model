# MAGIC HACK: without this, a recent Homebrew update broke GDAL!
GDAL_LIBRARY_PATH = "/usr/local/lib/libgdal.dylib"
import ctypes
ctypes.CDLL(GDAL_LIBRARY_PATH)

from rasterstats import zonal_stats
import glob
import os
import csv
import datetime
import fiona
import multiprocessing
import pandas

def parse_filename(filename):
    fn = os.path.split(filename)[1]
    year = fn[3:7]
    month = fn[7:9]
    return year, month

def parse_buffer_radius(filename):
    name = os.path.split(filename)[1]
    return name[-7:9]


class RunningAvg:
    def __init__(self):
        self.count = 0
        self.average = 0

    def update(self, value):
        self.count += 1
        if self.count == 1:
            self.average = value
        else:
            self.average = ((self.count/(self.count+1)) * self.average) + (1/(self.count+1) * value)


class PB:
    """
    Phils useful percent done tracker class!
    Initiate with total iterations
    .pd() returns percent done and updates the count
    .eta() returns an estimaed time until done
    """
    def __init__(self,length):
        self.i = 0
        self.length = length
        self.t1 = datetime.datetime.now()
        self.td = RunningAvg()

    def pd(self,inc=1):
        p = 100 * self.i/self.length
        self.i += inc
        self.td.update(datetime.datetime.now() - self.t1)
        self.t1 = datetime.datetime.now()
        return str(p)[0:4]

    def eta(self):
        return str(datetime.datetime.now() + (self.length - self.i)*self.td.average)

def extract_worker(raster):
    year, month = parse_filename(raster)
    stats = zonal_stats(shapefile, raster, stats=stat)
    return [[site, statistic[stat], month, year] for site, statistic in zip(shapefile_order, stats)]



if __name__ == "__main__":
    # Read command line arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--buffers", help="Directory to buffer shapefiles")
    parser.add_argument("-r", "--rasters", help="path to rasters")
    parser.add_argument("-o", "--output", help="directory for output csv file")
    parser.add_argument("-s", "--stat", help="statistic of ['mean','min','max','median','majority','minority']")
    parser.add_argument("-j", "--jobs")
    args = parser.parse_args()
    shapefile_path = args.buffers or "Buffers/"
    output_path = args.output or "SiteLevelSIC/"
    stat = args.stat or 'mean'
    raster_path = args.rasters or "PassiveMicrowaveSIC/"
    jobs = args.jobs or 1
    jobs = int(jobs)
    rasters = glob.glob(raster_path + '*.tif')
    if os.path.isdir(shapefile_path):
        buffer_files = glob.glob(shapefile_path + '*.shp')
    elif os.path.isfile(shapefile_path):
        buffer_files = [shapefile_path]
    else:
        raise IOError("Buffer files not found | {}".format(shapefile_path))

    print("Processing {} shapefiles and {} rasters".format(len(buffer_files), len(rasters)))
    # Get order of site_ids in shapefile
    for shapefile in buffer_files:
        r = parse_buffer_radius(shapefile)
        output = os.path.join(output_path,"ADPE_SIC_{}.csv".format(r))
        print("Processing radius {}km to {}".format(r,output))
        with fiona.open(shapefile) as sites:
            shapefile_order = [site['properties']['site_id'] for site in sites]


        # Multiprocessing using extract_worker()
        if jobs > 1:
            pool = multiprocessing.Pool(jobs)
            data = pool.map(extract_worker,rasters)
            with open('seaice_vals_tmp.csv','w') as csvfile:
                writer = csv.writer(csvfile)
                for raster in data:
                    for shape in raster:
                        writer.writerow(shape)
        else:
            # Progress status bar
            pb = PB(len(rasters))
            with open('seaice_vals_tmp.csv','a') as csvfile:
                writer = csv.writer(csvfile)

                # Iterate though input rasters and extract statistics for each polygon
                for raster in rasters:
                    year, month = parse_filename(raster)
                    stats = zonal_stats(shapefile, raster, stats=stat)

                    # Recombine stats with site_ids and write intermediary long form output
                    for site, statistic in zip(shapefile_order, stats):
                        writer.writerow([site, statistic[stat], month, year])
                    print("Procesed {}/{}  {}% complete  ETA {}".format(month, year, pb.pd(), pb.eta()))

        # Reload intermediary datafile
        data = pandas.read_csv('seaice_vals_tmp.csv',names=['site', 'SIC', 'month', 'year'])
        # Fix year to MAPPPD season/year form
        season_mod = [x <= 5 for x in data.month]
        data['season'] = pandas.Series(data.year - season_mod, index=data.index)

        # Rescale raster to percent
        data.SIC = (data.SIC / 250) * 100

        # Long to wide format
        data = data.pivot_table(index=['site', 'season'], values='SIC', columns='month')
        data = data.reset_index()
        data.columns = ['site_id', 'year', 'SIC_MONTH_1', 'SIC_MONTH_2', 'SIC_MONTH_3', 'SIC_MONTH_4',
                        'SIC_MONTH_5', 'SIC_MONTH_6', 'SIC_MONTH_7', 'SIC_MONTH_8', 'SIC_MONTH_9',
                        'SIC_MONTH_10', 'SIC_MONTH_11', 'SIC_MONTH_12']

        # Write output
        data.to_csv(output)
