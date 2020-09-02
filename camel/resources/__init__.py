from pathlib import Path

DIR_RESOURCES = Path(__file__).parent.absolute()

# CSS
CSS_STYLE = DIR_RESOURCES / 'css' / 'style.css'

# FONTS
FONT_SANS = DIR_RESOURCES / 'fonts' / 'FreeSans.ttf'
FONT_SANS_BOLD = DIR_RESOURCES / 'fonts' / 'FreeSansBold.ttf'

# YAML
YAML_ASSEMBLY_VELVETOPTIMISER = DIR_RESOURCES / 'yaml', 'assembly_velvetoptimiser.yml'
YAML_QUALITY_CHECKS = DIR_RESOURCES / 'yaml' / 'quality_checks.yml'
YAML_SAMTOOLS_VARIANT_CALLING = DIR_RESOURCES / 'yaml' / 'samtools_variant_calling.yml'
YAML_READ_MAPPING_BOWTIE2 = DIR_RESOURCES / 'yaml' / 'read_mapping_bowtie2.yml'
YAML_READ_TRIMMING = DIR_RESOURCES / 'yaml' / 'read_trimming.yml'
YAML_RES_CHAR_FAST = DIR_RESOURCES / 'yaml' / 'resistance_characterization_fast.yml'
YAML_TRIMMING_ASSEMBLY_SPADES = DIR_RESOURCES / 'yaml' / 'trimming_assembly_spades.yml'
YAML_TYPING_FAST = DIR_RESOURCES / 'yaml' / 'sequence_typing_fast.yml'

# Static
LOGO_SCIENSANO = DIR_RESOURCES / 'static' / 'logo-sciensano.png'

# Citations
DIR_CITATIONS = DIR_RESOURCES / 'citations'
