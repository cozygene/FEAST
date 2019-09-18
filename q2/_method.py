import os
import qiime2
import pandas as pd
import tempfile
import subprocess


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
    missing_ = list(set(fmeta[source_sink_column]) \
                    - set(source_ids + sink_ids))
    if len(missing_) > 0:
        raise ValueError(('The sub-classes of %s must'
                          ' be given as a source or sink: %s')\
                          %(str(source_sink_column),
                          ', '.join(map(str, missing_))))
    # rename ids in source and sink columns (only if all rep.)
    rename_ = {**{id_:'Source' for id_ in source_ids},
                **{id_:'Sink' for id_ in sink_ids}}
    fmeta[source_sink_column].replace(to_replace=rename_,
                                            inplace=True)

    return fmeta


def sourcetrack(table: pd.DataFrame,
                metadata: pd.DataFrame,
                environment_column: str,
                source_sink_column: str,
                source_ids: list,
                sink_ids: list,
                shared_id_column: str = None,
                EM_iterations: int = 1000,
                different_sources: bool = True) -> pd.DataFrame:

    # split the ids used for sources and sinks
    source_ids = source_ids.split(",")
    sink_ids = sink_ids.split(",")

    # create metadata formatted for FEAST
    # check if there are shared ids.
    if shared_id_column is not None:
        keep_cols = [environment_column,
                     source_sink_column,
                     sink_ids]
        rename_cols =  ['Env', 'SourceSink', 'id']
    else:
        keep_cols = [environment_column,
                     source_sink_column]
        rename_cols = ['Env', 'SourceSink']
    # keep only those columns
    feast_meta = metadata.dropna(subset=keep_cols)[keep_cols]
    # filter the metadata & table so they are matched
    shared_index = list(set(table.columns) & set(feast_meta.index))
    feast_meta = feast_meta.reindex(shared_index)
    table = table[shared_index]
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
    if different_sources == True:
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

        # build command for feast 
        cmd = ['FEAST.R',
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
