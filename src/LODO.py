import subprocess
base_command = [
    "python", "Main.py",
    "--raw-data", r"./BioMetaPipe/Result/data/16s/Raw/otu/Raw_feature_finally.csv",
    "--meta-data", r"./BioMetaPipe/Result/data/16s/Raw_meta_finally.csv",
    "--feature-type", "16s",
    "--groups", "CTR", "ADA",
    "--class-level", "otu",
    "--filter-threshold","0.0001",
    "--mapping", "CTR=0", "ADA=1",
    "--filter-mode", "abundance",
    "--output-path", r"./BioMetaPipe/Result/",
    "--feature","wilcoxon",
    "--norm-mode","std",
    "--correction", "no_correct","--no-feature-selection",
    "--no-cv"
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