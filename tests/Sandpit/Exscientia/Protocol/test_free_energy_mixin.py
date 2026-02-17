import numpy as np
import pandas as pd
import pytest

from BioSimSpace.Sandpit.Exscientia.Protocol import _FreeEnergyMixin


class Test_List_num_lam:
    """Test if the new _FreeEnergyMixin class could handle the num_lam as input
    with input as float.
    """

    @staticmethod
    @pytest.fixture(scope="class")
    def protocol():
        protocol = _FreeEnergyMixin(
            lam=0.5, lam_vals=None, min_lam=0.0, max_lam=1.0, num_lam=11
        )
        return protocol

    def test_getLambdaValues(self, protocol):
        value = protocol.getLambdaValues(type="list")
        assert isinstance(value, list)
        assert isinstance(value[5], float)
        assert np.isclose(value[5], 0.5)

    def test_getLambda(self, protocol):
        value = protocol.getLambda(type="float")
        assert isinstance(value, float)
        assert np.isclose(value, 0.5)

    def test_getLambdaIndex(self, protocol):
        value = protocol.getLambdaIndex()
        assert isinstance(value, int)
        assert np.isclose(value, 5)

    def test_getLambdaValues_df(self, protocol):
        value = protocol.getLambdaValues(type="dataframe")
        assert isinstance(value, pd.DataFrame)
        assert isinstance(value.loc[5], pd.Series)
        assert np.isclose(value.loc[5], 0.5)

    def test_getLambda_Series(self, protocol):
        value = protocol.getLambda(type="series")
        assert isinstance(value, pd.Series)
        assert np.isclose(value, 0.5)


class Test_List_lam_vals(Test_List_num_lam):
    """Test if the new _FreeEnergyMixin class could handle the lam_vals as input
    with input as float.
    """

    @staticmethod
    @pytest.fixture(scope="class")
    def protocol():
        protocol = _FreeEnergyMixin(
            lam=0.5,
            # Note that 0.5 is the first to check if being sorted
            lam_vals=[0.5, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0],
        )
        return protocol


class Test_df_lam_vals(Test_List_lam_vals):
    """Test if the new _FreeEnergyMixin class could handle the lam_vals as input
    with input as pd.DataFrame.
    """

    @staticmethod
    @pytest.fixture(scope="class")
    def protocol():
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"fep": 0.5}),
            # Note that 0.5 is the first to check if being sorted
            lam_vals=pd.DataFrame(
                data={"fep": [0.5, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]}
            ),
        )
        return protocol


class Test_df_num_lam(Test_List_num_lam):
    """Test if the new _FreeEnergyMixin class could handle the num_lam as input
    with input as pd.DataFrame.
    """

    @staticmethod
    @pytest.fixture(scope="class")
    def protocol():
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"fep": 0.5}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": 0.0}),
            max_lam=pd.Series(data={"fep": 1.0}),
            num_lam=11,
        )
        return protocol


def test_df_num_lam():
    """Test using multiple column df as input to num_lam."""
    protocol = _FreeEnergyMixin(
        lam=pd.Series(data={"bonded": 0.5, "coul": 0.5, "vdw": 0.5}),
        lam_vals=None,
        min_lam=pd.Series(data={"bonded": 0.0, "coul": 0.0, "vdw": 0.0}),
        max_lam=pd.Series(data={"bonded": 1.0, "coul": 1.0, "vdw": 1.0}),
        num_lam=11,
    )
    assert protocol.getLambdaIndex() == 5


def test_df_num_lam_different_index():
    """Test the case where min_lam and max_lam don't have the same index."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"bonded": 0.5}),
            lam_vals=None,
            min_lam=pd.Series(data={"bonded": 0.0}),
            max_lam=pd.Series(data={"coul": 1.0}),
            num_lam=11,
        )


def test_larger_max_lam():
    """Test the case where max_lam larger than 1."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"fep": 0.5}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": 0.0}),
            max_lam=pd.Series(data={"fep": 1.1}),
            num_lam=11,
        )


def test_smaller_min_lam():
    """Test the case where min_lam smaller than 0."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"fep": 0.5}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": -0.1}),
            max_lam=pd.Series(data={"fep": 1.0}),
            num_lam=11,
        )


def test_smaller_num_lam():
    """Test the case where num_lam smaller than 2."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"fep": 0.5}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": 0.1}),
            max_lam=pd.Series(data={"fep": 1.0}),
            num_lam=1,
        )


def test_invert_min_max():
    """Test the case where min_lam is larger than max_lam."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"fep": 0.5}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": 1.0}),
            max_lam=pd.Series(data={"fep": 0.0}),
            num_lam=11,
        )


class Test_df_lam_vals:
    """Test the case where a 'real' multiple column lambda dataframe is given."""

    @staticmethod
    @pytest.fixture(scope="class")
    def protocol():
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"bonded": 0.75, "coul": 0.5, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "bonded": [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                    "coul": [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
        )
        return protocol

    def test_getLambdaIndex(self, protocol):
        assert protocol.getLambdaIndex() == 3

    def test_getLambda(self, protocol):
        """Check that a warning is given when the Lambda (with multiple fields)
        cannot be converted to a float.
        """
        with pytest.warns(UserWarning):
            Lambda = protocol.getLambda(type="float")
            assert isinstance(Lambda, pd.Series)

    def test_getLambdaValues(self, protocol):
        """Check that a warning is given when the LambdaValues (with multiple
        columns) cannot be converted to a list of float.
        """
        with pytest.warns(UserWarning):
            Lambda = protocol.getLambdaValues(type="list")
            assert isinstance(Lambda, pd.DataFrame)


def test_not_in_num_lam():
    """Test The case where the lambda is not in lambda values generated by num_lam."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"bonded": 0.55, "coul": 0.55, "vdw": 0.55}),
            lam_vals=None,
            min_lam=pd.Series(data={"bonded": 0.0, "coul": 0.0, "vdw": 0.0}),
            max_lam=pd.Series(data={"bonded": 1.0, "coul": 1.0, "vdw": 1.0}),
            num_lam=11,
        )


def test_num_lam_not_int():
    """Test The case where the num_lam is not int."""
    with pytest.raises(TypeError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(
                data={
                    "bonded": 0.55,
                }
            ),
            lam_vals=None,
            min_lam=pd.Series(
                data={
                    "bonded": 0.0,
                }
            ),
            max_lam=pd.Series(
                data={
                    "bonded": 1.0,
                }
            ),
            num_lam=0.1,
        )


def test_not_in_lam_vals():
    """Test The case where the lambda is not in lambda values."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"bonded": 0.6, "coul": 0.5, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "bonded": [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                    "coul": [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
        )


def test_not_in_range():
    """Test The case where the lambda values is not in the range [0,1]."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(
                data={
                    "bonded": 0.6,
                }
            ),
            lam_vals=pd.DataFrame(
                data={
                    "bonded": [-0.1, 0.6],
                }
            ),
        )


def test_duplicated():
    """Test The case where lambda values has duplication."""
    with pytest.raises(ValueError):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"bonded": 0.75, "coul": 0.5, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "bonded": [0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                    "coul": [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
        )


def test_nonvalid_index_name():
    """Test The case where lambda name aaa is not in the permitted list."""
    with pytest.warns(UserWarning):
        protocol = _FreeEnergyMixin(
            lam=pd.Series(data={"aaa": 0.75, "coul": 0.5, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "aaa": [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                    "coul": [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
        )
