# Configuration file where you can set the parameter default values and
# descriptions.
DEFAULT_SHARED = None
DEFAULT_EMITR = 1000
DEFAULT_DIFFS = True

DESC_META = ('Sample metadata file containing sources'
             ' and sinks for source tracking.')
DESC_TBL = ('Feature table file containing sources'
            ' and sinks for source tracking.')
DESC_MP = ('The mixing proportions returned from FEAST.'
           ' The mixing proportions table is an S1 by '
           'S2 matrix P, where S1 is the number sinks '
           'and S2 is the number of sources (including '
           'an unknown source). Each row in matrix P sums'
           ' to 1. Pij is the contribution of source j to '
           'sink i. If Pij == NA it indicates that source '
           'j was not used in the analysis of sink i.')
DESC_ENVC = ('Sample metadata column with a description '
             'of the sampled environment (e.g., human gut).')
DESC_SSC = ('Sample metadata column with labels for source or a sink.'
            'All the sub-classes in this column must be in'
            ' either source_ids or sink_ids.')
DESC_SOURCEID = (
    'Comma-separated list (without spaces) of class ids '
    'contained in source_sink_column to be considered as sources.')
DESC_SINKID = ('Comma-separated list (without spaces) of class ids '
               'contained in source_sink_column to be considered as sinks.')
DESC_SHARED = ('Sample metadata column with the Sink-Source id.'
               ' When using multiple sinks, each tested with the '
               'same group of sources')
DESC_EMITR = ('A numeric value indicating the number of EM iterations.')
DESC_DIFFS = ('A Boolean value indicating the source-sink assignment.'
              'Different-sources is True if different sources are assigned'
              'to each sink, otherwise different-sources should be False.')
