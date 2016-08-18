import os
from os.path import abspath, dirname, join, isdir
from distutils.dir_util import copy_tree
import pandas as pd
import biom
from biom.util import compute_counts_per_sample_stats
import seaborn as sns
import matplotlib.pyplot as plt
import statistics

SKIN = "green"

def otu_table(output_dir: str, table: biom.Table) -> None:

    # Set output directory location
    output = join(output_dir, 'q2-summarizer-resources/')

    # Get summary stats
    results, df = _get_summary_stats(table)

    # Set range as a varibale
    range_value = '{}-{}'.format(results['minimum'], results['maximum'])

    # Move HTML resources to directory
    current_file_path = abspath(__file__)
    current_dir_path = dirname(current_file_path)
    copy_tree(join(current_dir_path, "dist/"), join(output, "dist/"))

    # Make histogram plot and save to disk
    sns.set_style("whitegrid")
    sns_plot = sns.distplot(results['depths'], kde = False,
                                               rug = True,
                                               hist = True,
                                               color='#3498db')
    sns_plot.set(xlabel = 'Sampling depths', ylabel = 'Frequency')
    fig = sns_plot.get_figure()
    fig.savefig(join(output, "histogram.png"), dpi = 300)
    fig.savefig(join(output, "histogram.pdf"), dpi = 300)

    # Write HTML file to disk
    with open(join(output_dir, "index.html"), 'w') as html:
        html.write(_header)
        html.write(_new_row())
        html.write(_make_sbox(title= 'Number of samples', value = results['samples'], color = "red"))
        html.write(_make_sbox(title= 'Number of OTUs', value = results['total_otu'], color = "orange"))
        html.write(_make_sbox(title= 'Total sequencing depth', value = results['total_counts'], color = "yellow"))
        html.write(_make_sbox(title= 'Denisty of non-zeros', value = results['density'], color = "green"))
        html.write(_end())
        html.write(_new_row())
        html.write(_make_sbox(title = 'Mean', value = results['mean']))
        html.write(_make_sbox(title = 'Median', value = results['median']))
        html.write(_make_sbox(title= 'Range', value = range_value))
        html.write(_make_sbox(title= 'Standard deviation', value = results['std']))
        html.write(_end())
        html.write(_plots)
        html.write(_make_datatable(df))
        html.write(_footer)


def _get_summary_stats(table):

    minimum, maximum, median, mean, counts_per_sample = \
        compute_counts_per_sample_stats(table)

    depths = list(counts_per_sample.values())
    std = round(statistics.stdev(depths), 2)
    total_otu = int(sum(depths))
    samples = int(len(depths))
    total_counts = int(table.length(axis="observation"))
    density = round(table.get_table_density(), 2) * 100
    depths = list(counts_per_sample.values())

    depths_table = pd.Series(counts_per_sample, name='Depths')
    depths_table.index.name = 'SampleID'
    depths_table = depths_table.reset_index()


    # Place results in a dictionary
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

    # Return data needed for plots
    return(results, depths_table)


def _make_sbox(title = 'Title', value = 18, icon = 'ion-pie-graph', color = 'blue'):
    return('''
            <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-''' + color + '''">
                    <div class="inner">
                        <h3>''' + str(value) + '''</h3>
                        <p>''' + str(title) + '''</p>
                    </div>
                    <div class="icon">
                        <i class="ion ''' + icon + '''"></i>
                    </div>
                </div>
            </div>
    ''')

def _make_datatable(table):
    return('''
            <div class="row">
              <div class="col-md-12">
                <div class="box box-solid">
                  <div class="box-header with-border">
                    <h3 class="box-title">Sampling Depths Table</h3>
                  </div>
                  <div class="box-body">
                    <div class="dataTables_wrapper form-inline dt-bootstrap">
                     ''' + table.to_html(index = False, classes = ['table', 'table-bordered', 'table-hover', 'datatables']) + '''
                    </div>
                  </div>
                </div>
              </div>
            </div>
    ''')

def _new_row():
    return('''<div class="row">''')

def _end():
    return('''</div>''')


_header = '''
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <title>QIIME2-Summariser</title>
  <meta content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no" name="viewport">
  <link rel="stylesheet" href="q2-summarizer-resources/dist/bootstrap/css/bootstrap.min.css">
  <link rel="stylesheet" href="q2-summarizer-resources/dist/plugins/datatables/jquery.dataTables.min.css">
  <link rel="stylesheet" href="q2-summarizer-resources/dist/css/font-awesome.min.css">
  <link rel="stylesheet" href="q2-summarizer-resources/dist/css/ionicons.min.css">
  <link rel="stylesheet" href="q2-summarizer-resources/dist/css/AdminLTE.min.css">
  <link rel="stylesheet" href="q2-summarizer-resources/dist/css/skins/_all-skins.min.css">
</head>
<body class="hold-transition skin-''' + SKIN + ''' layout-top-nav">
<div class="wrapper">
  <header class="main-header">
    <nav class="navbar navbar-static-top">
      <div class="container">
        <div class="navbar-header">
          <a href="http://qiime.org/" class="navbar-brand">
            <b>QIIME2-</b>Summarizer
          </a>
          <button type="button"
                  class="navbar-toggle collapsed"
                  data-toggle="collapse"
                  data-target="#navbar-collapse">
            <i class="fa fa-bars"></i>
          </button>
        </div>
        <div class="collapse navbar-collapse pull-left" id="navbar-collapse">
          <ul class="nav navbar-nav">
            <li>
                <a href="#preprocessing">FASTQ data</a>
            </li>
            <li>
                <a href="#preprocessing">Joined PE reads</a>
            </li>
            <li>
                <a href="#preprocessing">Split Libraries (de-multiplexing)</a>
            </li>
            <li>
                <a href="#preprocessing">OTU-Picking</a>
            </li>
            <li class="active">
            <a href="#otu">OTU Table <span class="sr-only">(current)</span></a>
            </li>
          </ul>
        </div>
      </div>
    </nav>
  </header>
  <div class="content-wrapper">
    <div class="container">
        <section class="content">
'''

_plots = '''
<div class="row">
  <div class="col-md-6">
    <div class="box box-solid">
      <div class="box-header with-border">
        <h3 class="box-title">Sampling Depth Histogram</h3>
        <a href="q2-summarizer-resources/histogram.png" class="btn btn-primary btn-sm pull-right" type="button">PNG</a>
        <a href="q2-summarizer-resources/histogram.pdf" class="btn btn-primary btn-sm pull-right" type="button">PDF</a>
      </div>
      <div class="box-body">
         <img class="img-responsive pad" src="q2-summarizer-resources/histogram.png">
      </div>
    </div>
  </div>
  <div class="col-md-6">
    <div class="box box-solid">
      <div class="box-header with-border">
        <h3 class="box-title">OTU Rank Abundance</h3>
      </div>
      <div class="box-body">
         <img class="img-responsive pad" src="q2-summarizer-resources/rank_abundance.png">
      </div>
    </div>
  </div>
</div>
'''

_footer = '''
</section>
</div>
</div>
<footer class="main-footer">
<div class="container">
<div class="pull-right hidden-xs">
  <strong>
    <a href="http://almsaeedstudio.com">Made with Admin LTE 2.3.6</a>
  </strong>
</div>
</div>
</footer>
</div>

<script src="q2-summarizer-resources/dist/plugins/jQuery/jquery-2.2.3.min.js"></script>
<script src="q2-summarizer-resources/dist/plugins/datatables/dataTables.bootstrap.min.js"></script>
<script src="q2-summarizer-resources/dist/plugins/datatables/jquery.dataTables.min.js"></script>
<script src="q2-summarizer-resources/dist/bootstrap/js/bootstrap.min.js"></script>
<script src="q2-summarizer-resources/dist/js/app.min.js"></script>

<script>
$(function () {
    $(".datatables").DataTable();});
</script>
</body>
</html>
'''
