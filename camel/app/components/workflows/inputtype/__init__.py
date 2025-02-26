from camel.app.components.workflows.inputtype.fastahelper import FastaHelper
from camel.app.components.workflows.inputtype.illuminahelper import IlluminaHelper
from camel.app.components.workflows.inputtype.onthelper import ONTHelper

helper_by_input_type = {
    'fasta': FastaHelper,
    'illumina': IlluminaHelper,
    'ont': ONTHelper
}
