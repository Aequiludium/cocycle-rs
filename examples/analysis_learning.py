"""Validate and classify the fixed graph CSV from diagram_analysis --export."""

import argparse
import csv
import math
from pathlib import Path
import sys

COLUMNS = ["sample_id", "split", "label", "source", "vertex_count", "edge_count",
           "mean_degree", "leaf_fraction"]


def read_samples(path):
    """Preserve row identity/order and reject incomplete or mismatched fixtures."""
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != COLUMNS:
            raise ValueError("unexpected CSV columns or feature order")
        rows = list(reader)
    expected_ids = {f"{label}-{n}" for n in (3, 4, 6, 5, 7)
                    for label in ("path", "cycle")}
    seen = set()
    for row in rows:
        sample_id = row["sample_id"]
        if sample_id not in expected_ids or sample_id in seen:
            raise ValueError("unknown or duplicate sample_id")
        seen.add(sample_id)
        label, size = sample_id.split("-")
        n = int(size)
        split = "train" if n in (3, 4, 6) else "test"
        if (row["label"] != label or row["split"] != split
                or row["source"] != "unweighted-graph-v1"):
            raise ValueError("sample identity, label, split or source mismatch")
        values = [float(row[column]) for column in COLUMNS[4:]]
        # Independent combinatorial expectations, not another exporter call.
        expected = ([n, n - 1, 2 * (n - 1) / n, 2 / n] if label == "path"
                    else [n, n, 2, 0])
        if any(not math.isfinite(value) or abs(value - wanted) > 1e-14
               for value, wanted in zip(values, expected)):
            raise ValueError("nonfinite or inconsistent graph features")
        if None in row:
            raise ValueError("extra CSV values")
    if seen != expected_ids:
        raise ValueError("empty or incomplete sample set")
    return rows


def learn(rows):
    # Optional dependencies are loaded only by this explicitly requested example.
    import numpy as np
    import sklearn
    import scipy
    import joblib
    import threadpoolctl
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.tree import DecisionTreeClassifier

    train = [row for row in rows if row["split"] == "train"]
    test = [row for row in rows if row["split"] == "test"]
    # IDs, labels, splits and counts never enter the feature matrix.
    x_train = np.array([[float(row["leaf_fraction"])] for row in train])
    x_test = np.array([[float(row["leaf_fraction"])] for row in test])
    model = Pipeline([("scale", StandardScaler()),
                      ("tree", DecisionTreeClassifier(max_depth=1, random_state=0))])
    model.fit(x_train, [row["label"] for row in train])
    scaler = model.named_steps["scale"]
    if scaler.n_samples_seen_ != 6 or abs(scaler.mean_[0] - 0.25) > 1e-14:
        raise ValueError("preprocessing did not use only the declared training rows")
    predicted = model.predict(x_test).tolist()
    wanted = [row["label"] for row in test]
    if predicted != wanted:
        raise ValueError(f"teaching predictions differ: {predicted} != {wanted}")
    print(f"scikit-learn={sklearn.__version__} numpy={np.__version__} "
          f"scipy={scipy.__version__} joblib={joblib.__version__} "
          f"threadpoolctl={threadpoolctl.__version__}")
    print("train=6 test=4; training-only scaler mean=0.25")
    for row, prediction in zip(test, predicted):
        print(f"{row['sample_id']}: expected={row['label']} predicted={prediction}")
    print("4/4 fixed teaching labels matched; no application accuracy claim.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    args = parser.parse_args(argv)
    try:
        learn(read_samples(args.csv))
    except (OSError, ValueError, TypeError, ImportError) as error:
        print(f"Learning example failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
