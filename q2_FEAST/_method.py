import os
import tempfile
import subprocess
import pandas as pd
from ._feast_defaults import (DEFAULT_DIFFS,
                              DEFAULT_EMITR)


def run_commands(cmds, verbose=True):
    """
    This function is a script runner.
    It was obtained from https://github.com/ggloor
    /q2-aldex2/blob/master/q2_aldex2/_method.py
    """
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def feast_format(fmeta: pd.DataFrame,
                 source_sink_column: str,
                 source_ids: list,
                 sink_ids: list,) -> pd.DataFrame:
    """
    Helper function to format metadata for FEAST.
    """

    # ensure that all sub-cats in SourceSink are represented
    missing_ = list(set(fmeta[source_sink_column])
                    - set(source_ids + sink_ids))
    if len(missing_) > 0:
        raise ValueError(('All of the sub-classes of %s must'
                          ' be given as a source or sink'
                          'the sub-class(es) [%s] are missing.')
                         % (str(source_sink_column),
                             ', '.join(map(str, missing_))))
    # rename ids in source and sink columns (only if all rep.)
    rename_ = {**{id_: 'Source' for id_ in source_ids},
               **{id_: 'Sink' for id_ in sink_ids}}
    fmeta[source_sink_column].replace(to_replace=rename_,
                                      inplace=True)

    return fmeta


def sourcetrack(table: pd.DataFrame,
                metadata: pd.DataFrame,
                environment_column: str,
                source_sink_column: str,
                source_ids: list,
                sink_ids: list,
                shared_id_column: str,
                em_iterations: int = DEFAULT_EMITR,
                different_sources: bool = DEFAULT_DIFFS) -> pd.DataFrame:

    # split the ids used for sources and sinks
    source_ids = source_ids.split(",")
    sink_ids = sink_ids.split(",")

    # create metadata formatted for FEAST
    # check if there are shared ids.
    # currently FEAST requires an id
    # column but in future versions it will
    # be an optional peram.
    if shared_id_column is not None:
        keep_cols = [environment_column,
                     source_sink_column,
                     shared_id_column]
        rename_cols = ['Env', 'SourceSink', 'id']
    else:
        keep_cols = [environment_column,
                     source_sink_column]
        rename_cols = ['Env', 'SourceSink']

    # import and check all columns given are in dataframe
    metadata = metadata.to_dataframe()
    metadata.index = metadata.index.astype(str)
    if not all([col_ in metadata.columns for col_ in keep_cols]):
        raise ValueError('Not all columns given are present in the'
                         ' sample metadata file. Please check that'
                         ' the input columns are in the given metdata.')

    # keep only those columns
    feast_meta = metadata.dropna(subset=keep_cols)
    feast_meta = feast_meta.loc[:, keep_cols]

    # filter the metadata & table so they are matched
    table = table.T
    shared_index = list(set(table.columns) & set(feast_meta.index))
    feast_meta = feast_meta.reindex(shared_index)
    table = table.loc[:, shared_index]

    # format the sub-classes for source-sink
    feast_meta = feast_format(feast_meta,
                              source_sink_column,
                              source_ids,
                              sink_ids)
    if shared_id_column is not None:
        # encode the shared SourceSink id column
        #  with numerics ranging from 1-N
        shared_ = set(metadata[shared_id_column])
        rename_ = {id_: i for i, id_ in enumerate(shared_)}
        feast_meta[shared_id_column].replace(to_replace=rename_,
                                             inplace=True)
    # rename those columns for FEAST
    feast_meta.columns = rename_cols

    # if there are different sources
    if different_sources:
        different_sources = 1
    else:
        different_sources = 0

    # save all intermediate files into tmp dir
    with tempfile.TemporaryDirectory() as temp_dir_name:
        # save the tmp dir locations
        biom_fp = os.path.join(temp_dir_name, 'input.tsv.biom')
        map_fp = os.path.join(temp_dir_name, 'input.map.txt')
        summary_fp = os.path.join(temp_dir_name, 'output.proportions.txt')

        # Need to manually specify header=True for Series (i.e. "meta"). It's
        # already the default for DataFrames (i.e. "table"), but we manually
        # specify it here anyway to alleviate any potential confusion.
        table.to_csv(biom_fp, sep='\t', header=True)
        feast_meta.to_csv(map_fp, sep='\t', header=True)

        # build command for FEAST
        cmd = ['source_tracking.R',
               biom_fp,
               map_fp,
               different_sources,
               summary_fp]
        cmd = list(map(str, cmd))

        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            raise Exception("An error was encountered while running FEAST"
                            " in R (return code %d), please inspect stdout"
                            " and stderr to learn more." % e.returncode)

        # if run was sucessfull import the data and return
        proportions = pd.read_csv(summary_fp, index_col=0)
        proportions.index.name = "sampleid"
        return proportions
