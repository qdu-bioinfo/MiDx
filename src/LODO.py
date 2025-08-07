import subprocess
base_command = [
    "python", "./src/Main.py",
    "--raw-data", r"./data/wgs/Raw/t_sgb/Raw_feature_finally.csv",
    "--meta-data", r"./data/wgs/Raw/Raw_meta_finally.csv",
    "--feature-type", "wgs",
    "--groups", "CTR", "CRC",
    "--class-level", "t_sgb",
    "--filter-threshold","0.0001",
    "--mapping", "CTR=0", "CRC=1",
    "--filter-mode", "abundance",
    "--output-path", r"./Result/",
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