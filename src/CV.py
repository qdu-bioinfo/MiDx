import subprocess
import os
base_command = [
    "python", "./src/Main.py",
    # Input raw feature table (CSV), samples × features
    "--raw-data", r"./data/wgs/Raw/t_sgb/Raw_feature_finally.csv",
    # Input sample metadata (CSV), index is sample IDs
    "--meta-data", r"./data/wgs/Raw/Raw_meta_finally.csv",
    # Feature type tag (e.g., 16s, wgs)
    "--feature-type", "wgs",
    # Groups to include (here: control CTR and adenoma ADA)
    "--groups", "CTR", "CRC",
    # Classification level tag (using OTU level here)
    "--class-level", "t_sgb",
    # Filter threshold: retain features with abundance >= 0.0001
    "--filter-threshold", "0.0001",
    # Label mapping: map CTR to 0, ADA to 1
    "--mapping", "CTR=0", "CRC=1",
    # Filter mode: based on abundance
    "--filter-mode", "abundance",
    # Output directory for results
    "--output-path", r"./Result/",
    # Feature selection method: Wilcoxon test
    "--feature", "wilcoxon",
    # Normalization mode: standardization (std)
    "--norm-mode", "std",
    # Batch correction method: no correction ("no_correct" is interpreted as the correction parameter value)
    "--correction", "no_correct",
    # Skip the feature selection stage
    "--no-feature-selection",
    # Skip the LODO model stage
    "--no-lodo"
]
cmd = base_command.copy()
result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace"
            )
print("Exit code:", result.returncode)
if result.stdout:
    print("Output:\n", result.stdout)
if result.stderr:
    print("Errors:\n", result.stderr)

print("\nAll commands executed!")