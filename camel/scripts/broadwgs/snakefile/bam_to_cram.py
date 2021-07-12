from pathlib import Path

SNAKEFILE_BAMTOCRAM = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_bamtocram = Path('bamtocram')
OUTPUT_BAMTOCRAM_CRAM = _dir_bamtocram / "convert" / "cram.io"
OUTPUT_BAMTOCRAM_CRAM_checksum = _dir_bamtocram / "checksum" / 'cram.md5'
OUTPUT_BAMTOCRAM_CRAI = _dir_bamtocram / "index" /  'crai.io'
OUTPUT_BAMTOCRAM_CRAM_metrics = _dir_bamtocram / "metrics" / 'cram_validation_report.io'