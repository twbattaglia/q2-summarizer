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
    name='Visualize and summarize data',
    description='Visualize data from an OTU table and generate plots related to'
                'sampling depths and OTU abundances.'
)
