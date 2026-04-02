"""
Moran's Eigenvector Maps (MEM) Generator for Spatial Autocorrelation Analysis

Implements the PCNM (Principal Coordinates of Neighbour Matrices) methodology
from Borcard & Legendre (2002) to generate spatial eigenvectors that can be
used as predictors in regression/gradient forest models to account for
spatial autocorrelation.

References:
    - Borcard, D., & Legendre, P. (2002). All-scale spatial analysis of ecological
      data by means of principal coordinates of neighbour matrices. Ecological
      Modelling, 153(1-2), 51-68.
    - Fitzpatrick, M. C., & Keller, S. R. (2015). Ecological genomics meets
      community-level modelling of biodiversity: mapping the genomic landscape
      of current and future environmental adaptation. Ecology Letters, 18(1), 1-16.

Usage:
    python scripts/generate_mems.py --input data.csv --output data_with_mems.csv \\
        --lat decimalLatitude --lon decimalLongitude

    # Include all positive eigenvalue MEMs
    python scripts/generate_mems.py --input data.csv --output out.csv --all

    # Specify exact number of MEMs
    python scripts/generate_mems.py --input data.csv --output out.csv --n-mems 5
"""

import argparse
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix


def calculate_morans_i(y: np.ndarray, coords: np.ndarray) -> float:
    """
    Calculate Moran's I statistic for spatial autocorrelation.

    Moran's I ranges from -1 (perfect dispersion) to +1 (perfect clustering),
    with 0 indicating random spatial distribution.

    Args:
        y: Array of values to test for spatial autocorrelation (n_samples,)
        coords: Array of coordinates (n_samples, 2) as [lat, lon]

    Returns:
        Moran's I statistic value
    """
    y = np.array(y)
    coords = np.array(coords)
    n = len(y)

    # Compute distance matrix and inverse-distance weights
    dists = distance_matrix(coords, coords)

    with np.errstate(divide='ignore'):
        weights = 1.0 / dists
    np.fill_diagonal(weights, 0)

    # Row-normalize weights
    row_sums = weights.sum(axis=1)
    row_sums[row_sums == 0] = 1
    weights = weights / row_sums[:, np.newaxis]

    # Calculate Moran's I
    y_mean = np.mean(y)
    y_diff = y - y_mean

    denominator = np.sum(y_diff ** 2)
    if denominator == 0:
        return 0.0

    numerator = np.sum(weights * np.outer(y_diff, y_diff))
    morans_i = (n / np.sum(weights)) * (numerator / denominator)

    return morans_i


def generate_mems(coords: np.ndarray, truncation_factor: float = 4.0) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Generate Moran's Eigenvector Maps (MEMs) / PCNM variables.

    Implements the Borcard & Legendre (2002) methodology:
    1. Compute Euclidean distance matrix D
    2. Find truncation threshold via Minimum Spanning Tree (max MST edge)
    3. Construct truncated spatial weighting matrix W
    4. Gower-center the matrix
    5. Eigendecomposition
    6. Return eigenvectors with positive eigenvalues

    Args:
        coords: Array of coordinates (n_samples, 2) as [lat, lon]
        truncation_factor: Multiplier for threshold in weighting function (default 4.0)
                          Values >= 4 yield identical results up to multiplicative constant

    Returns:
        Tuple of:
            - DataFrame with MEM_1, MEM_2, ... columns (all positive eigenvalue eigenvectors)
            - Array of corresponding eigenvalues (sorted descending)
    """
    coords = np.array(coords)
    n = len(coords)

    # Step 1: Euclidean distance matrix
    D = distance_matrix(coords, coords)

    # Step 2: Minimum Spanning Tree to find truncation threshold
    # MST ensures all points remain connected at the threshold distance
    D_sparse = csr_matrix(np.triu(D))
    mst = minimum_spanning_tree(D_sparse)
    mst_array = mst.toarray()
    threshold = np.max(mst_array)

    print(f"Number of samples: {n}")
    print(f"Truncation threshold (max MST edge): {threshold:.4f}")

    # Step 3: Construct truncated spatial weighting matrix
    # Weighting function from Borcard & Legendre (2002): w_ij = 1 - (d_ij / 4t)^2
    W = np.zeros((n, n))
    mask = (D <= threshold) & (D > 0)  # Neighbors within threshold, exclude self
    W[mask] = 1 - (D[mask] / (truncation_factor * threshold)) ** 2

    # Step 4: Gower centering: W_centered = H @ W @ H where H = I - 11'/n
    H = np.eye(n) - np.ones((n, n)) / n
    W_centered = H @ W @ H

    # Step 5: Eigendecomposition (eigh for symmetric matrices)
    eigvals, eigvecs = np.linalg.eigh(W_centered)

    # Sort by eigenvalues in descending order
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Step 6: Select positive eigenvalues only
    # Negative eigenvalues correspond to complex/non-Euclidean components
    pos_indices = eigvals > 1e-9  # Numeric tolerance for zero
    mem_vectors = eigvecs[:, pos_indices]
    mem_values = eigvals[pos_indices]

    print(f"Generated {len(mem_values)} MEMs with positive eigenvalues")

    # Create DataFrame with MEM columns
    col_names = [f"MEM_{i+1}" for i in range(len(mem_values))]
    mem_df = pd.DataFrame(mem_vectors, columns=col_names)

    return mem_df, mem_values


def add_mems_to_data(
    input_csv: str,
    output_csv: str,
    lat_col: str = "decimalLatitude",
    lon_col: str = "decimalLongitude",
    n_mems: int | None = None,
    use_all: bool = False
) -> pd.DataFrame:
    """
    Load data, generate MEMs, and save combined output.

    By default, uses the first half of MEMs with positive eigenvalues
    as recommended by Fitzpatrick & Keller (2015).

    Args:
        input_csv: Path to input CSV file
        output_csv: Path to output CSV file
        lat_col: Name of latitude column
        lon_col: Name of longitude column
        n_mems: Explicit number of MEMs to include (overrides default)
        use_all: If True, include all positive eigenvalue MEMs

    Returns:
        DataFrame with original data plus MEM columns
    """
    # Load data
    print(f"Loading data from: {input_csv}")
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} samples with {len(df.columns)} columns")

    # Validate coordinate columns exist
    if lat_col not in df.columns:
        raise ValueError(f"Latitude column '{lat_col}' not found in data")
    if lon_col not in df.columns:
        raise ValueError(f"Longitude column '{lon_col}' not found in data")

    # Extract coordinates
    coords = df[[lat_col, lon_col]].values

    # Check for missing coordinates
    valid_mask = ~np.isnan(coords).any(axis=1)
    if not valid_mask.all():
        n_invalid = (~valid_mask).sum()
        print(f"Warning: {n_invalid} samples have missing coordinates and will have NaN MEMs")

    # Generate MEMs using valid coordinates
    valid_coords = coords[valid_mask]
    mem_df, mem_values = generate_mems(valid_coords)

    # Determine number of MEMs to use
    total_mems = len(mem_values)
    if use_all:
        n_to_use = total_mems
    elif n_mems is not None:
        n_to_use = min(n_mems, total_mems)
    else:
        # Default: first half of positive eigenvalue MEMs (Fitzpatrick & Keller 2015)
        n_to_use = total_mems // 2
        if n_to_use < 1 and total_mems > 0:
            n_to_use = 1

    print(f"Using {n_to_use} of {total_mems} positive eigenvalue MEMs")

    # Subset MEMs
    mem_df = mem_df.iloc[:, :n_to_use]

    # Handle samples with missing coordinates
    if not valid_mask.all():
        # Create full MEM dataframe with NaNs for invalid rows
        full_mem_df = pd.DataFrame(
            np.nan,
            index=range(len(df)),
            columns=mem_df.columns
        )
        full_mem_df.loc[valid_mask, :] = mem_df.values
        mem_df = full_mem_df

    # Reset index to align with original dataframe
    mem_df.index = df.index

    # Combine original data with MEMs
    final_df = pd.concat([df, mem_df], axis=1)

    # Save output
    final_df.to_csv(output_csv, index=False)
    print(f"Saved output to: {output_csv}")
    print(f"Output has {len(final_df.columns)} columns ({n_to_use} MEM columns added)")

    return final_df


def main():
    parser = argparse.ArgumentParser(
        description="Generate Moran's Eigenvector Maps (MEMs) for spatial autocorrelation analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Default: use first half of positive eigenvalue MEMs
    python generate_mems.py --input data.csv --output data_with_mems.csv

    # Include all positive eigenvalue MEMs
    python generate_mems.py --input data.csv --output data_with_mems.csv --all

    # Specify exact number of MEMs to include
    python generate_mems.py --input data.csv --output data_with_mems.csv --n-mems 5

References:
    Borcard & Legendre (2002). Ecological Modelling 153:51-68
    Fitzpatrick & Keller (2015). Ecology Letters 18:1-16
        """
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input CSV file path"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output CSV file path"
    )
    parser.add_argument(
        "--lat",
        default="decimalLatitude",
        help="Name of latitude column (default: decimalLatitude)"
    )
    parser.add_argument(
        "--lon",
        default="decimalLongitude",
        help="Name of longitude column (default: decimalLongitude)"
    )
    parser.add_argument(
        "--n-mems",
        type=int,
        default=None,
        help="Number of MEMs to include (default: first half of positive eigenvalue MEMs)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        dest="use_all",
        help="Include all MEMs with positive eigenvalues"
    )

    args = parser.parse_args()

    add_mems_to_data(
        input_csv=args.input,
        output_csv=args.output,
        lat_col=args.lat,
        lon_col=args.lon,
        n_mems=args.n_mems,
        use_all=args.use_all
    )


if __name__ == "__main__":
    main()
