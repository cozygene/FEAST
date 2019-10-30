import os
import unittest
import pandas as pd
from biom import Table
from qiime2 import Artifact, Metadata
from q2_FEAST._method import microbialtracking
from numpy.testing import assert_allclose


class TestFEAST(unittest.TestCase):

    def setUp(self):

        in_dir = os.path.dirname(os.path.abspath(__file__))
        # path goes to tutorial data
        def get_file(dir_, file_): return os.path.join(dir_,
                                                       os.pardir,
                                                       'tutorials',
                                                       'data',
                                                       'backhed',
                                                       file_)
        # path goes to tutorial data
        def get_test(dir_, file_): return os.path.join(dir_,
                                                       'data',
                                                       file_)

        # import multi test data
        q2bt_multi = Artifact.load(get_file(in_dir, 'table-multi.qza'))
        self.q2bt_multi = q2bt_multi.view(Table).to_dataframe().T
        self.q2mf_multi = Metadata.load(
            get_file(in_dir, 'metadata-multi.qza'))
        # perams
        self.envcol = 'Env'
        self.sourcesinkcol = 'SourceSink'
        self.sourceids = 'Source'
        self.sinkids = 'Sink'
        self.sharedidcolumn = 'host_subject_id'
        # expected
        self.exp_mixed = pd.read_csv(get_test(in_dir, 'exp-mixed-mp.tsv'),
                                     sep='\t',
                                     index_col=0)

    def test_sourcetrack_multi(self):

        # run FEAST mixed source tracker (False is slow)
        res_mpdf = microbialtracking(self.q2bt_multi,
                                     self.q2mf_multi,
                                     self.envcol,
                                     self.sourcesinkcol,
                                     self.sourceids,
                                     self.sinkids,
                                     self.sharedidcolumn,
                                     different_sources = True,
                                     em_iterations = 5000)
        #tres_all.to_csv('/Users/cameronmartino/Downloads/testing-feast.tsv', sep='\t')                           
        res_ = res_mpdf.loc['Unknown',
                            self.exp_mixed.columns].values
        exp_ = self.exp_mixed.loc['Unknown', :].values
        assert_allclose(res_, exp_, atol=.50)

    def test_sourcetrack_value_errors(self):

        # test missing id
        with self.assertRaises(ValueError):
            microbialtracking(self.q2bt_multi,
                              self.q2mf_multi,
                              self.envcol,
                              self.sourcesinkcol,
                              'Sour',
                              self.sinkids,
                              self.sharedidcolumn)

        # test bad column
        with self.assertRaises(ValueError):
            microbialtracking(self.q2bt_multi,
                              self.q2mf_multi,
                              'oh-no-bad-bad',
                              self.sourcesinkcol,
                              self.sourceids,
                              self.sinkids,
                              self.sharedidcolumn,
                              em_iterations=100)


if __name__ == "__main__":
    unittest.main()
