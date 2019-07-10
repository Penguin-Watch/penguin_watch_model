#MODIFIED FROM MAPPPD REPO - written by Chris Che-Castaldo (https://github.com/CCheCastaldo)

# MAGIC HACK: without this, a recent Homebrew update broke GDAL!
#GDAL_LIBRARY_PATH = "/usr/local/lib/libgdal.dylib"
#import ctypes
#ctypes.CDLL(GDAL_LIBRARY_PATH)

import fiona
from fiona import collection
from fiona.crs import from_epsg
from shapely.geometry import mapping, shape, MultiPolygon

def largest(polygon):
    '''
    return index of largest multipolgyon
    :param polygon: shapely.geometry.MultiPolygon
    :return: index: int
    '''
    return max(range(len(polygon)), key = lambda x: polygon[x].area)

def buffer_validity_test(buffers):
    '''
    checks polygons for non-null geometry.
    returns tuple of boolean validity and invalid polygons (or empty list)
    :param polygons:
    :return:
    '''
    validity = [x.type != 'Polygon' for x in buffers]
    return not any(validity),[i for i,x in enumerate(validity)if x]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sites", help="shapefile of colony locations (EPSG:3031)")
    parser.add_argument("-r", "--radius", nargs='+', help="radius of buffers")
    parser.add_argument("-o", "--output", help="output shapefile")
    parser.add_argument("-l", "--land", help="shapefile of land mask")
    args = parser.parse_args()

    site_locations = args.sites or 'Locations/SitesEPSG3031.shp'
    land_mask = args.land or 'LandMask/land.shp'
    try:
        radius_list = [int(args.radius)]
    except TypeError:
        radius_list = [int(x) for x in args.radius]
    output_path = args.output or 'Buffers/'
    with fiona.open(land_mask) as land_mask_shape:
        land = MultiPolygon([shape(lm['geometry']) for lm in land_mask_shape])
        land = land.buffer(0)
    for radius in radius_list:

        output = output_path + 'buffer{0:03d}.shp'.format(radius)
        print("Processing radius {}km to {}".format(radius,output))

        props = []
        points = []

        with collection(site_locations, "r") as input:
            schema = input.schema.copy()
            for point in input:
                points.append(shape(point['geometry']))
                props.append(point['properties'])
            buffers = [point.buffer(radius*1000).difference(land) for point in points]
            for ind, buffer in enumerate(buffers):
                if buffer.type == 'MultiPolygon':
                    buffers[ind] = buffer[largest(buffer)]
            schema = {'geometry': 'Polygon', 'properties': {'site_id': 'str'}}
            validity,issues = buffer_validity_test(buffers)
            if not validity:
                print("Buffer @ {}km skipped due to NULL geometry. Try a larger radius.".format(radius))
                print("No data found in range of {}".format([props[i]['site_id'] for i in issues]))
            else:
                with collection(
                        output, "w", "ESRI Shapefile", schema,crs=from_epsg(3031)) as output:
                    for ind, buffer in enumerate(buffers):
                        output.write({
                            'properties': {
                                'site_id': props[ind]['site_id']
                            },
                            'geometry': mapping(buffer)
                        })
