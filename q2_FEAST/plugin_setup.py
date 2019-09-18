
import qiime2.sdk
import qiime2.plugin
from qiime2.plugin import (Properties, Int,
                           Metadata, Str, Bool)
from q2_types.feature_table import (FeatureTable,
                                    Frequency, 
                                    RelativeFrequency)

from ._method import sourcetrack
from ._feast_defaults import (DESC_META, DESC_ENVC,
                              DESC_SSC, DESC_SOURCEID,
                              DESC_SINKID, DESC_SHARED,
                              DESC_EMITR, DESC_DIFFS,
                              DESC_MP, DESC_TBL)

# TODO: will need to fix the version number
__version__ = '0.1.0'

# perams types
PARAMETERS = {'metadata': Metadata,
              'environment_column': Str,
              'source_sink_column': Str,
              'source_ids': Str,
              'sink_ids': Str,
              'shared_id_column': Str,
              'EM_iterations': Int,
              'different_sources': Bool}
# perams descriptions
PARAMETERDESC = {'metadata': DESC_META,
                 'environment_column': DESC_ENVC,
                 'source_sink_column': DESC_SSC,
                 'source_ids': DESC_SOURCEID,
                 'sink_ids': DESC_SINKID,
                 'shared_id_column': DESC_SHARED,
                 'EM_iterations': DESC_EMITR,
                 'different_sources': DESC_DIFFS}

citations = qiime2.plugin.Citations.load('citations.bib',
                                         package='q2_FEAST')

plugin = qiime2.plugin.Plugin(
    name='FEAST',
    version=__version__,
    website="https://github.com/cozygene/FEAST",
    citations=[citations['Shenhav2019-ca']],
    short_description=('Plugin for FEAST source-tracking'),
    description=('This is a QIIME 2 plugin supporting microbial'
                 ' source-tracking through FEAST.'),
    package='q2_FEAST')

plugin.methods.register_function(
    function=sourcetrack,
    inputs={'table': FeatureTable[Frequency]},
    parameters=PARAMETERS,
    outputs=[('mixing_proportions', FeatureTable[RelativeFrequency])],
    input_descriptions={'table': DESC_TBL},
    parameter_descriptions=PARAMETERDESC,
    output_descriptions={'mixing_proportions': DESC_MP},
    name='microbial source-tracking',
    description=('A major challenge of analyzing the compositional '
                'structure of microbiome data is identifying its potential origins. '
                'Here, we introduce fast expectation-maximization microbial source '
                'tracking (FEAST), a ready-to-use scalable framework that can '
                'simultaneously estimate the contribution of thousands of potential '
                'source environments in a timely manner, thereby helping unravel the '
                'origins of complex microbial communities. The information gained from '
                'FEAST may provide insight into quantifying contamination, tracking the '
                'formation of developing microbial communities, as well as distinguishing '
                'and characterizing bacteria-related health conditions.'),
)
