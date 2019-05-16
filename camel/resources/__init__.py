import os
_resources_folder = os.path.dirname(os.path.realpath(__file__))

# CSS
CSS_STYLE = os.path.join(_resources_folder, 'css', 'style.css')

# FONTS
FONT_SANS = os.path.join(_resources_folder, 'fonts', 'FreeSans.ttf')
FONT_SANS_BOLD = os.path.join(_resources_folder, 'fonts', 'FreeSansBold.ttf')

# YAML
YAML_ASSEMBLY_VELVETOPTIMISER = os.path.join(_resources_folder, 'yaml', 'assembly_velvetoptimiser.yml')
YAML_QUALITY_CHECKS = os.path.join(_resources_folder, 'yaml', 'quality_checks.yml')
YAML_SAMTOOLS_VARIANT_CALLING = os.path.join(_resources_folder, 'yaml', 'samtools_variant_calling.yml')
YAML_READ_MAPPING_BOWTIE2 = os.path.join(_resources_folder, 'yaml', 'read_mapping_bowtie2.yml')
YAML_READ_TRIMMING = os.path.join(_resources_folder, 'yaml', 'read_trimming.yml')
YAML_RES_CHAR_FAST = os.path.join(_resources_folder, 'yaml', 'resistance_characterization_fast.yml')
YAML_TRIMMING_ASSEMBLY_SPADES = os.path.join(_resources_folder, 'yaml', 'trimming_assembly_spades.yml')
YAML_TYPING_FAST = os.path.join(_resources_folder, 'yaml', 'sequence_typing_fast.yml')

# Static
LOGO_SCIENSANO = os.path.join(_resources_folder, 'static', 'logo-sciensano.png')
