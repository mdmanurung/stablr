# Ledoit-Wolf and equicorrelated S-matrix fixture

These frozen values were generated with Python 3.13.3, NumPy 2.2.4,
scikit-learn 1.7.1, and knockpy 1.3.5. `LedoitWolf` used its default
`assume_centered=False` and `block_size=1000` settings. The S-matrix used
`knockpy.smatrix.compute_smatrix(..., method="equicorrelated")` with its
default tolerance `1e-5` and ten-iteration PSD scaling search.

The fixture tests numerical semantics only. It is independent of the OOL
release packet and cannot be used as promotion evidence.
