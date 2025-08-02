import subprocess
base_command = [
    "python", "Main.py",
    "--raw-data", r"./BioMetaPipe/Result/data/wgs/Raw/t_sgb/Raw_feature_finally.csv",
    "--meta-data", r"./BioMetaPipe/Result/data/wgs/Raw_meta_finally.csv",
    "--feature-type", "wgs",
    "--groups", "CTR", "ADA",
    "--class-level", "t_sgb",
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