from camel.app.components.workflows.readtype.illuminahelper import IlluminaHelper
from camel.app.components.workflows.readtype.onthelper import ONTHelper

helper_by_input_type = {
    'illumina': IlluminaHelper,
    'ont': ONTHelper
}
