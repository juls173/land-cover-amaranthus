import pandas as pd
import numpy as np

# Read the CSV file
df = pd.read_csv("agri_gf_input.csv")

# Identify SNP columns (those starting with "AgCMH_")
snp_columns = [col for col in df.columns if col.startswith("AgCMH_")]

print(f"Found {len(snp_columns)} SNP columns to impute")
print(f"Total rows: {len(df)}")

# Track columns to drop
columns_to_drop = []

# Impute missing values for each SNP column with the column mean
for col in snp_columns:
    if df[col].isnull().any():
        count = df[col].isnull().sum()
        missing_proportion = count / len(df[col])
        
        # Drop columns with MORE than 10% missing (>0.1, not <0.1)
        if missing_proportion > 0.1:
            columns_to_drop.append(col)
            print(f"Column {col} marked for dropping ({missing_proportion:.1%} missing)")
            continue
        
        # Calculate mean of non-missing values
        col_mean = df[col].mean()
        
        # Impute with mean (fixed syntax)
        df[col].fillna(col_mean, inplace=True)
        
        print(f"{col}: Imputed {count} missing values with mean {col_mean:.4f}")

# Drop all marked columns at once
if columns_to_drop:
    df.drop(columns=columns_to_drop, inplace=True)
    print(f"\nDropped {len(columns_to_drop)} columns with >10% missing data")

# Update snp_columns list to reflect dropped columns
snp_columns = [col for col in snp_columns if col not in columns_to_drop]

# Save the imputed data
df.to_csv('agri_imputed_data.csv', index=False)

print("\nImputation complete! Saved to 'agri_imputed_data.csv'")

# Summary statistics
print("\nSummary:")
print(f"Total SNP columns processed: {len(snp_columns)}")
print(f"Data shape: {df.shape}")
print(f"Remaining missing values in SNP columns: {df[snp_columns].isnull().sum().sum()}")