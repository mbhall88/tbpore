from pathlib import Path

TMP_NAME = ".tbpore"
repo_root = Path(__file__).parent.parent.resolve()
H37RV_genome = repo_root / "data/H37RV_genome/h37rv.fa.gz"
H37RV_mask = repo_root / "data/H37RV_genome/compass-mask.bed"
cache_dir = repo_root / ".cache"
config_file = repo_root / ".config.yaml"
external_scripts_dir = repo_root / "external_scripts"
