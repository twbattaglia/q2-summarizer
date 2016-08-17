import qiime.plugin
import q2_summarizer
from ._otu_table import otu_table
from q2_types import (FeatureTable, Frequency)

plugin = qiime.plugin.Plugin(
    name='summarizer',
    version=q2_summarizer.__version__,
    website='https://github.com/twbattaglia/q2-summarizer',
    package='q2_summarizer',
    user_support_text=None,
    citation_text=None
)

plugin.visualizers.register_function(
    function=otu_table,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    name='Summarize OTU-table data.',
    description='This visualizer produces an HTML visualization of two '
                'key-value mappings, each sorted in alphabetical order by key.'
)
