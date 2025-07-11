import logging
import time

import ee


def get_info(ee_obj, max_retries=4):
    """Make an exponential back off getInfo call on an Earth Engine object"""
    # output = ee_obj.getInfo()
    output = None
    for i in range(1, max_retries):
        try:
            output = ee_obj.getInfo()
        except ee.ee_exception.EEException as e:
            if ('Earth Engine memory capacity exceeded' in str(e) or
                    'Earth Engine capacity exceeded' in str(e) or
                    'Too many concurrent aggregations' in str(e) or
                    'Computation timed out.' in str(e)):
                # TODO: Maybe add 'Connection reset by peer'
                logging.info(f'    Resending query ({i}/{max_retries})')
                logging.info(f'    {e}')
            else:
                # TODO: What should happen for unexpected EE exceptions?
                #   It might be better to reraise the exception and exit
                logging.info(f'    {e}')
                logging.info('    Unhandled Earth Engine exception')
                continue
        except Exception as e:
            logging.info(f'    Resending query ({i}/{max_retries})')
            logging.debug(f'    {e}')

        if output is not None:
            break

        time.sleep(i ** 3)

    return output


def constant_image_value(image, crs='EPSG:32613', scale=1):
    """Extract the output value from a "constant" image"""
    rr_params = {
        'reducer': ee.Reducer.first(),
        'geometry': ee.Geometry.Rectangle([0, 0, 10, 10], crs, False),
        'scale': scale,
    }
    return get_info(ee.Image(image).reduceRegion(**rr_params))


def point_image_value(image, xy, scale=1):
    """Extract the output value from a calculation at a point"""
    rr_params = {
        'reducer': ee.Reducer.first(),
        'geometry': ee.Geometry.Point(xy),
        'scale': scale,
    }
    return get_info(ee.Image(image).reduceRegion(**rr_params))
