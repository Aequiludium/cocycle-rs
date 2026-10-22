# Persistent spectral features and external learning

[Documentation](../index.md) / Research

This page describes the reproduction under review in
[PR #129](https://github.com/Aequiludium/cocycle-rs/pull/129), checked against
[source 71aa864](https://github.com/Aequiludium/cocycle-rs/blob/71aa864cd3b7afafe541ba5e88c5866c90d175a8/tools/reproduce_persistent_spectral_features.py).
It is absent from this documentation branch. It is a standalone Python
reproduction using existing Rust descriptor workers, with no new public Rust
spectral-feature or classifier interface. Production integration with the
[persistent Laplacian candidate](../guides/persistent-laplacian.md) remains separate.

## From pair spectra to a feature row

For a finite increasing grid $t_0,\ldots,t_{m-1}$, the reproduction constructs
real unweighted persistent Laplacians for pairs $i\le j$ in degrees zero and
one. It computes complete tiny spectra independently in Python, selects the
smallest eigenvalues and retains multiplicities. Ordinary Hodge spectra and
Rust Betti/lifetime descriptors are separate controls.

`vectorize_spectra(samples, grid, degrees=(0, 1), length=8,
operator="persistent_laplacian")` returns owned arrays and metadata. This is a
function in the reproduction script, not a crate API.

| Read | Contract |
| --- | --- |
| `values` | Real array of shape `(sample_count, degree_count * m*(m+1)/2 * length)` |
| `validity_mask` | Boolean array of the same shape; true for computed entries, false for exhausted-domain padding |
| `columns` | Degree, source-grid index, target-grid index, scales and eigenvalue rank, in that order |
| `sample_ids`, `source_ids`, `bases` | Unique sample identities, declared sources and ordered bases for interpretation |
| `domain_sizes`, `spectrum_coverage`, `numeric_contract` | Per-pair domain/coverage and shared units, real inner product and numerical origin |

With two degrees, four grid points and length eight, each row has 160 columns.
Padding is permitted only when the domain has fewer than `length` eigenvalues;
failed or missing solves are errors. Input records declare complete spectra or
an adequate smallest prefix. The function validates that declaration's shape and
provenance, but does not independently certify the producer's extremal indexes.
The residual-based partial Hodge/Dirac results explicitly leave index coverage
uncertified and cannot simply be relabeled as a certified smallest prefix.

Values within the recorded numerical zero tolerance are set to zero before
sorting. A valid zero may therefore be exact or thresholded; the mask distinguishes
computed values from padding, not exact zeros from tiny eigenvalues. Truncation
and padding do not normalize complex size. A truncated prefix sum is not the trace.

## Reproduce the experiment

Run these commands in the candidate checkout with a new output directory.
The author file must be the unchanged `src/perslap.py` from
[b36b080](https://github.com/tomogwen/persistentlaplaciandatascience/blob/b36b080f739653d84284f3f7ca9ccea5fad286ae/src/perslap.py);
the script checks the pinned source identity.

```sh
python -m pip install numpy==2.2.6 scipy==1.15.3 gudhi==3.12.0 scikit-learn==1.7.2 threadpoolctl==3.6.0
python -m unittest discover -s tools -p test_reproduce_persistent_spectral_features.py
python tools/reproduce_persistent_spectral_features.py --author-source PATH_TO_PERSLAP --output target/spectral-feature-replay --seed 20261012
```

Python dependencies belong to this external workflow. Rust descriptor workers
compile and run against the candidate's existing library; they do not introduce
Python or learning dependencies into the crate. The output directory is never
overwritten. It retains source/environment hashes, generated fixtures, schemas,
spectra, splits, predictions, validation outcomes and failure records outside Git.

## Interpret the controls and evidence

Two group-held-out synthetic experiments each have 64 matched pairs, with 96
training rows and 32 test rows. Betti/lifetime and ordinary spectral controls
score 50% in both. Persistent degree zero separates the tree controls at 100%,
where shortest-path/effective-resistance controls already separate them.
Persistent degree one scores 93.75% on the disk controls, whose different cycle
fillings have spectra $\{3,3,3\}$ and $\{1,3,3\}$ despite matching PH and degree-zero
spectra. The forests use values only; mask-only, size and values-plus-mask
ablations remain separate work. Shared motifs occur in both splits, so these
numbers do not establish new-family or dataset-level generalization.

Core checks compare 160 worked pair matrices using restricted-boundary SVD,
exact rational ranks, pinned author algebra and declared GUDHI topology cases.
The original feature loop runs with audited callbacks. This finite evidence does
not identify finite-field topology with real rank on arbitrary complexes.
The inclusive uniform grid differs from the author's uneven grid; the original
MNIST/MoleculeNet pipelines and degree-two experiments are not reproduced here.

The [resource study](https://github.com/Aequiludium/cocycle-rs/blob/71aa864cd3b7afafe541ba5e88c5866c90d175a8/benches/reports/spectral-feature-resources.md)
retains ten fresh complete workflows and a discarded warmup. Its wall time
includes tracemalloc, fresh Rust child builds, checks and both forests; parent
Python working-set peak is not whole-process-tree memory or algorithm allocation.
These are resource snapshots, not spectral-kernel rankings or production learning
acceptance. Further descriptor/stability work must preserve source, scale,
dimension, spectral coverage and the assumptions of the cited theorem.
