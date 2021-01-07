import pytest
import numpy as np
import pandas as pd

@pytest.fixture
def sample_maf():
    df = pd.read_csv("examples/data/sample_from_maf.tsv", sep="\t")
    df = df.replace(["NA", "na", "-", "Unknown"], np.nan)
    return df

