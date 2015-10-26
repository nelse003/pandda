
import numpy

from scitbx.array_family import flex

def calculate_occupancy_subtracted_map(ref_map, query_map, subtract_occ):
    """Calculate query_map - subtract_occ*ref_map"""
    assert isinstance(ref_map, flex.double), 'ref_map must of type flex.double'
    assert isinstance(query_map, flex.double), 'query_map must of type flex.double'
    assert ref_map.all() == query_map.all(), 'ref_map and query_map are not the same size: {!s} != {!s}'.format(ref_map.all(), query_map.all())
    assert ref_map.size() == query_map.size(), 'ref_map and query_map are not the same size: {!s} != {!s}'.format(ref_map.size(), query_map.size())
    # Important for subtraction to be "ref_map*subtract_occ" so that the output is a flex.double array
    diff_map = query_map - ref_map*subtract_occ
    return diff_map

def calculate_occupancy_correlations(ref_map, query_map, feature_region, min_occ=0.0, max_occ=1.0, occ_increment=0.01, reference_region=None, verbose=False):
    """Estimate the occupancy of a feature in "query", defined by grid points "feature_region". Reference should be an equivalent map differing only by not having the feature"""

    # Check that the inputs are the right types
    assert isinstance(ref_map, flex.double), 'ref_map must of type flex.double'
    assert isinstance(query_map, flex.double), 'query_map must of type flex.double'
    assert isinstance(feature_region, list), 'grid_point_selection must be of type list'
    assert isinstance(reference_region, list), 'grid_point_selection must be of type list'
    # Check that the two maps are the same size
    assert ref_map.all() == query_map.all(), 'ref_map and query_map are not the same size: {!s} != {!s}'.format(ref_map.all(), query_map.all())
    assert ref_map.size() == query_map.size(), 'ref_map and query_map are not the same size: {!s} != {!s}'.format(ref_map.size(), query_map.size())

    if verbose: print 'FEATURE REGION:', len(feature_region)
    if verbose: print 'REFERENCE REGION:', len(reference_region)

    ##########################################################################################

    # Calculate the map indices of the grid points where the feature is
    grid_accessor = ref_map.accessor()
    feature_region_idx = flex.size_t(map(grid_accessor, feature_region))
    # Calculate the map indices of the grid points of the reference region of the map
    if reference_region: reference_region_idx = flex.size_t(map(grid_accessor, reference_region))
    else:                reference_region_idx = None

    ##########################################################################################

    # Extract the reference map values around the site of interest
    feature_region_vals_ref = ref_map.as_1d().select(feature_region_idx)
    # Extract the reference map values for the correlation correction
    if reference_region: reference_region_vals_ref = ref_map.as_1d().select(reference_region_idx)
    else:                reference_region_vals_ref = ref_map.as_1d()

    ##########################################################################################

    return_vals = []

    # Create list of occuapncies to subtract
    all_feature_occs = numpy.arange(min_occ, max_occ, occ_increment)
    # Ensure last value is 1.0
    if all_feature_occs[-1] != max_occ:
        all_feature_occs = numpy.arange(min_occ, max_occ+occ_increment, occ_increment)
        all_feature_occs[-1] = max_occ

    # Iterate through different amounts of reference map subtraction to estimate feature occupancy
    for i_occ, feature_occ in enumerate(all_feature_occs):

        # Calculate how much of the reference map to take away
        ref_occ = 1 - feature_occ

        # Calculate the occupancy-subtracted difference map
        diff_map = calculate_occupancy_subtracted_map(ref_map=ref_map, query_map=query_map, subtract_occ=ref_occ)

        # Extract feature region values
        feature_region_vals_diff = diff_map.as_1d().select(feature_region_idx)
        # Extract reference region values
        if reference_region: reference_region_vals_diff = diff_map.as_1d().select(reference_region_idx)
        else:                reference_region_vals_diff = diff_map.as_1d()

        feature_region_corr     = flex.linear_correlation(feature_region_vals_diff, feature_region_vals_ref, subtract_mean=False)
        reference_region_corr   = flex.linear_correlation(reference_region_vals_diff, reference_region_vals_ref, subtract_mean=False)

        return_vals.append((feature_occ, feature_region_corr.coefficient(), reference_region_corr.coefficient()))

    return zip(*return_vals)

def calculate_maximum_series_discrepancy(labels, series_1, series_2):
    """Calculate the point at which two series are maximally different"""
    assert len(series_1) == len(series_2)
    assert len(series_1) == len(labels)
    diffs = [series_1[i] - series_2[i] for i in range(len(series_1))]
    return labels[diffs.index(max(diffs))]

def estimate_event_occupancy(ref_map, query_map, feature_region, min_occ=0.0, max_occ=1.0, occ_increment=0.01, reference_region=None, method='gradient', verbose=False):
    """Calculate an estimate for the occupancy of a defined event"""
    # Check method is valid
    assert method in ['gradient', 'value']
    # Calculate the correlations to the reference map for a range of occupancies
    event_occs, event_corrs, global_corrs = calculate_occupancy_correlations(
        ref_map          = ref_map,
        query_map        = query_map,
        feature_region   = feature_region,
        reference_region = reference_region,
        min_occ          = min_occ,
        max_occ          = max_occ,
        occ_increment    = occ_increment,
        verbose          = verbose)

    # Estimate the event occupancies
    if method == 'value':
        # Calculate the maximum discrepancy in the correlation values as an estimation of event occupancy
        occ_est = calculate_maximum_series_discrepancy(
            labels   = event_occs,
            series_1 = global_corrs,
            series_2 = event_corrs)
    elif method == 'gradient':
        # Calculate the maximum discrepancy in the correlation gradients as an estimation of event occupancy
        occ_est = calculate_maximum_series_discrepancy(
            labels   = event_occs[1:-1],
            series_1 = [event_corrs[i+1] - event_corrs[i-1] for i in range(1, len(event_corrs)-1)],
            series_2 = [global_corrs[i+1] - global_corrs[i-1] for i in range(1, len(global_corrs)-1)])

    return occ_est

