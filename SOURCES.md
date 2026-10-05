# Sources and reproducibility

The equations and numerical source locations are cited in the model and examples.
This public snapshot does not redistribute article PDFs, extracted article text,
copied figure images, unrelated MATLAB code or the full research archive.
Download papers from their publishers or open preprints when inspecting source
equations. Save the required PDFs under `sources/papers/` for historical source
validation scripts; their expected filenames are listed below.

Normal CLI use, the quick/full **unit** suites and the fast/reference grid
calibration run without PDFs. `scripts/validate/validate_stage3_a.py --full`
also hashes local article copies, so it requires the original paper files;
it is distinct from the self-contained full unit suite.

## Papers

- [bertaina2024.pdf](https://doi.org/10.1002/qute.202400032)
- [didomenico2010.pdf](https://doi.org/10.1364/AO.49.004801)
- [williams2008.pdf](https://doi.org/10.1364/JOSAB.25.001284)
- [kasdin1995.pdf](https://doi.org/10.1109/5.381848)
- [clivati2022.pdf](https://doi.org/10.1038/s41467-021-27808-1)
- [clivati2022-1.pdf](https://doi.org/10.1038/s41467-021-27808-1)
- [clivati2022-2.pdf](https://doi.org/10.1038/s41467-022-28447-w)
- [pittaluga2021.pdf](https://arxiv.org/abs/2012.15099)
- [zhou2023.pdf](https://arxiv.org/abs/2208.09347)
- [zeng2022_mode_pairing.pdf](https://doi.org/10.1038/s41467-022-31534-7)
- [zeng2022_mode_pairing_preprint.pdf](https://arxiv.org/abs/2201.04300)
- [zhu2023_mode_pairing.pdf](https://arxiv.org/abs/2208.05649)
- [zhou2023_async_mdi.pdf](https://arxiv.org/abs/2212.14190)
- [zhang2025_mode_pairing_prx.pdf](https://doi.org/10.1103/PhysRevX.15.021037)
- [xie2023_async_comparison.pdf](https://arxiv.org/abs/2302.14349)
- [jiang2008_urban.pdf](https://arxiv.org/abs/0807.1882)
- [li2013_lasers.pdf](https://thesis.caltech.edu/7799/)

## Reference software and data

[Bertaina et al., Zenodo record 10911121](https://doi.org/10.5281/zenodo.10911121)
provides `QKD.ipynb` and the four `*_meas.txt` files. The record declares
Creative Commons Attribution 4.0. Copies in `sources/data/` are unchanged.
See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Digitized curves and observed counters

`examples/b6/` and `sources/data/mp_phase/digitized/` contain numerical
digitizations, not copies of article figure images. Coordinates, units,
normalization decisions and source locations are retained in their JSON
manifests. They do not carry a new blanket MIT licence over the underlying
publications or third-party data. Mode-pairing exposure/count inputs cite the
source table for each experiment in `sources/data/mode_pairing/`.

## Evidence retained

Small Markdown/JSON records cover the independent implementation, Table I grid
convergence, the held-out TF comparison, mode-pairing phase analysis and the Zhu
diagnostics. Local absolute workspace paths have been made relative; numerical
results have not been changed. Large arrays, reports from each scan, logs and
the original Git history are deliberately omitted. No unpublished laboratory
observations are required by the supplied examples.
