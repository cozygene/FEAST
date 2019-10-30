import os
import shutil
import pandas as pd
from qiime2 import Artifact, Metadata
from q2_taxa._visualizer import barplot as _barplot

def barplot(output_dir: str,
            mixing_proportions: pd.DataFrame,
            metadata: Metadata) -> None:

    # fix identifiers
    mixing_proportions.index = [ind.split('_')[0]
                                for ind in mixing_proportions.index]
    mixing_proportions.columns = [col.split('_')[0]
                                  for col in mixing_proportions.columns]

    # import and check all columns given are in dataframe
    metadata = metadata.to_dataframe()

    # replace seperation character in metadata
    metadata = metadata.replace('_', '-',
                                regex=True)
    metadata.index = metadata.index.astype(str)
    metadata.index = [ind.replace('_', '-')
                    for ind in metadata.index]

    # make "Sink" metadata (columns of mixing_proportions)
    sink_metadata = metadata.reindex(mixing_proportions.columns)
    sink_metadata.index.name = 'sampleid'
    sink_metadata = Metadata(sink_metadata)

    # make "Source" metadata (index of mixing_proportions)
    source_metadata = metadata.reindex(mixing_proportions.index).fillna('Unkown-Source')
    source_metadata = source_metadata.astype(str).apply(lambda x: '; '.join(x), axis=1)
    source_metadata = pd.DataFrame(source_metadata,
                                columns = ['Taxon'])
    source_metadata.index.name = 'Feature ID'

    # make barplot
    _barplot(output_dir,
             mixing_proportions.T,
             pd.Series(source_metadata.Taxon),
             sink_metadata)

    # grab bundle location to fix
    bundle = os.path.join(output_dir,
                        'dist',
                        'bundle.js')
    # bundle terms to fix for our purpose
    bundle_rplc = {'Relative Frequency':'Source Contribution',
                'Taxonomic Level':'Source Grouping',
                'Sample':'Sink'}
    # make small text chnage to bundle
    with open(bundle) as f:
        newText=f.read()
        for prev, repl in bundle_rplc.items():
            newText = newText.replace(prev, repl)
    with open(bundle, "w") as f:
        f.write(newText)
