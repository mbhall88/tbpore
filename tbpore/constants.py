from pathlib import Path

TMP_NAME = ".tbpore"
repo_root = Path(__file__).parent.parent.resolve()
H37RV_genome = repo_root / "data/H37RV_genome/h37rv.fa.gz"
H37RV_mask = repo_root / "data/H37RV_genome/compass-mask.bed"
DECONTAMINATION_DB_INDEX = (
    Path.home() / ".tbpore/decontamination_db/remove_contam.map-ont.mmi"
)
DECONTAMINATION_DB_METADATA = repo_root / "data/decontamination_db/remove_contam.tsv.gz"
CACHE_DIR = Path.home() / ".cache"
CONFIG_FILE = repo_root / ".config.yaml"
EXTERNAL_SCRIPTS_DIR = repo_root / "external_scripts"
