from camel.app.components.workflows.readtype.illuminahelper import IlluminaHelper
from camel.app.components.workflows.readtype.iontorrenthelper import IonTorrentHelper
from camel.app.components.workflows.readtype.nanoporehelper import NanoporeHelper

helper_by_read_type = {
    'illumina': IlluminaHelper,
    'iontorrent': IonTorrentHelper,
    'nanopore': NanoporeHelper
}
