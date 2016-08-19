import os
from os.path import abspath, dirname, join, isdir
import jinja2
from jinja2 import Environment, FileSystemLoader, PackageLoader
from distutils.dir_util import copy_tree
import pandas as pd
import biom
from biom.util import compute_counts_per_sample_stats
import seaborn as sns
import matplotlib.pyplot as plt
import statistics


def otu_table(output_dir: str, table: biom.Table) -> None:

    output = join(output_dir, 'q2-summarizer-resources/')
    results, df = _get_summary_stats(table)
    range_value = '{}-{}'.format(results['minimum'], results['maximum'])

    # Move HTML resources to directory
    current_file_path = abspath(__file__)
    current_dir_path = dirname(current_file_path)
    copy_tree(join(current_dir_path, "dist/"), join(output, "dist/"))

    # Make depths dataframe html
    depth_df = df.to_html(index = False, classes = ['table', 'table-bordered', 'table-hover', 'datatables'])

    # Make histogram plot and save to disk
    sns.set_style("whitegrid")
    histogram_plot = sns.distplot(a = results['depths'],
                                  kde = False,
                                  rug = True,
                                  hist = True,
                                  color='#3498db')
    histogram_plot.set_title('Sampling Depths Distribution')
    histogram_plot.set_xlabel('Sampling Depth')
    histogram_plot.set_ylabel('Frequency')
    histogram_plot.get_figure().savefig(
        join(output, 'histogram.png'), dpi = 300)
    histogram_plot.get_figure().savefig(
        join(output, 'histogram.pdf'), dpi = 300)

    # Clear plots
    plt.gcf().clear()

    # OTU Rank abundance
    otuRA = ((table.sum(axis="observation") / results['total_counts'])/ 100)
    otu_rank = pd.DataFrame(data = otuRA, index = table.ids(axis='observation'))
    otu_rank.columns = ['Abundance']
    otu_rank_sorted = otu_rank.sort_values('Abundance', ascending = False).head(30)
    otu_rank_plot = sns.pointplot(x = otu_rank_sorted.index,
                                  y = otu_rank_sorted.Abundance,
                                  data = otu_rank_sorted)
    otu_rank_plot.set_title('OTU Rank Abundance')
    otu_rank_plot.set_yscale('log')
    otu_rank_plot.set_xlabel('Abundance Rank')
    otu_rank_plot.set_ylabel('Relative Abundance (log)')
    otu_rank_plot.set(xticklabels=[])
    otu_rank_plot.get_figure().savefig(
        join(output, 'rank_abundance.png'), dpi = 300)
    otu_rank_plot.get_figure().savefig(
        join(output, 'rank_abundance.pdf'), dpi = 300)

    # Add variables to jinja2 template
    env = Environment(loader = PackageLoader('q2_summarizer', 'templates'))
    template = env.get_template('otu-table.html')
    output_from_parsed_template = template.render(
        samples = str(results['samples']),
        total_otu = str(results['total_otu']),
        total_counts = str(results['total_counts']),
        density = str(results['density']),
        mean = str(results['mean']),
        median = str(results['median']),
        range_value = str(range_value),
        std = str(results['std']),
        table = depth_df
    )

    # Write HTML file to disk
    with open(join(output_dir, "index.html"), 'w') as html:
        html.write(output_from_parsed_template)

    return(None)


def _get_summary_stats(table):

    minimum, maximum, median, mean, counts_per_sample = \
        compute_counts_per_sample_stats(table)

    depths = list(counts_per_sample.values())
    std = round(statistics.stdev(depths), 2)
    total_counts = int(sum(depths))
    samples = int(len(depths))
    total_otu = int(table.length(axis="observation"))
    density = round(table.get_table_density(), 2) * 100
    depths = list(counts_per_sample.values())

    depths_table = pd.Series(counts_per_sample, name='Depths')
    depths_table.index.name = 'SampleID'
    depths_table = depths_table.reset_index()

    results = {'mean': int(mean),
               'median': int(median),
               'minimum': int(minimum),
               'maximum': int(maximum),
               'std': round(std, 2),
               'total_otu': int(total_otu),
               'samples': int(samples),
               'total_counts': int(total_counts),
               'density': str(density) + "%",
               'depths': depths}
    return(results, depths_table)
