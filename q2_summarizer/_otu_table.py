# This file is only included as an example. Remove this file when you are
# ready to develop your plugin.
import os
from os.path import abspath, dirname, join, isdir
from distutils.dir_util import copy_tree

import pandas as pd
import biom
from biom.util import compute_counts_per_sample_stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statistics


def otu_table(output_dir : str, table: biom.Table) -> None:

    # Set output directory location
    output = join(output_dir, 'q2-summarizer-resources/')
    rankabundance_output = join(output_dir, "rank_abundance.png")

    # Calculate all statistics
    minimum, maximum, median, mean, counts_per_sample = compute_counts_per_sample_stats(table)

    # Get sampling-depths values
    depths = list(counts_per_sample.values())

    # SD
    std = round(statistics.stdev(depths), 2)

    # Total OTU counts
    total_otu = int(sum(depths))
    samples = int(len(depths))

    # Number of different OTU's
    total_counts = int(table.length(axis="observation"))

    # Write HTML
    html = '''
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <meta http-equiv="X-UA-Compatible" content="IE=edge">
      <title>QIIME2-Summariser</title>
      <meta content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no" name="viewport">
      <link rel="stylesheet" href="q2-summarizer-resources/dist/bootstrap/css/bootstrap.min.css">
      <link rel="stylesheet" href="q2-summarizer-resources/dist/css/font-awesome.min.css">
      <link rel="stylesheet" href="q2-summarizer-resources/dist/css/ionicons.min.css">
      <link rel="stylesheet" href="q2-summarizer-resources/dist/css/AdminLTE.min.css">
      <link rel="stylesheet" href="q2-summarizer-resources/dist/css/skins/_all-skins.min.css">
    </head>
    <body class="hold-transition skin-blue layout-top-nav">
    <div class="wrapper">

      <header class="main-header">
        <nav class="navbar navbar-static-top">
          <div class="container">
            <div class="navbar-header">
              <a href="http://qiime.org/" class="navbar-brand"><b>QIIME2-</b>Summarizer</a>
              <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar-collapse">
                <i class="fa fa-bars"></i>
              </button>
            </div>

            <!-- Collect the nav links, forms, and other content for toggling -->
            <div class="collapse navbar-collapse pull-left" id="navbar-collapse">
              <ul class="nav navbar-nav">
                <li><a href="#preprocessing">FASTQ data</a></li>
                <li><a href="#preprocessing">Joined PE reads</a></li>
                <li><a href="#preprocessing">Split Libraries (de-multiplexing)</a></li>
                <li><a href="#preprocessing">OTU-Picking</a></li>
                <li class="active"><a href="#otu">OTU Table <span class="sr-only">(current)</span></a></li>
              </ul>
            </div>
          </div>
        </nav>

      </header>
      <!-- Full Width Column -->
      <div class="content-wrapper">
        <div class="container">

          <!-- Sampling depth -->
          <section class="content">
            <div class="row">
              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-red">
                  <div class="inner">
                    <h3>''' + str(round(samples, 0)) + '''</h3>
                    <p>Number of samples</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>

              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-orange">
                  <div class="inner">
                    <h3>''' + str(round(total_otu, 0)) + '''</h3>
                    <p>Number of OTU's</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>

              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-yellow">
                  <div class="inner">
                    <h3>''' + str(round(total_counts, 0)) + '''</h3>
                    <p>Total sequencing depth</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>

              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-green">
                  <div class="inner">
                    <h3>''' + 'coming soon' + '''</h3>
                    <p>Number of species</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>
            </div>

            <!-- Row 2 -->
            <div class="row">
              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-blue">
                  <div class="inner">
                    <h3>''' + str(round(int(mean), 0)) + '''</h3>
                    <p>Mean</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>

              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-blue">
                  <div class="inner">
                    <h3>''' + str(round(median, 0))  + '''</h3>
                    <p>Median</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>

              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-blue">
                  <div class="inner">
                    <h3>''' + str(int(minimum)) + '''-''' + str(int(maximum)) + '''</h3>
                    <p>Range</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>

              <div class="col-lg-3 col-xs-6">
                <div class="small-box bg-blue">
                  <div class="inner">
                    <h3>''' + str(std)  + '''</h3>
                    <p>Std Deviation</p>
                  </div>
                  <div class="icon">
                    <i class="ion ion-pie-graph"></i>
                  </div>
                </div>
              </div>
            </div>

            <!-- Rank abundance plot -->
            <div class="row">
              <div class="col-md-6">
                <div class="box box-solid">
                  <div class="box-header with-border">
                    <h3 class="box-title">Sampling Depth Histogram</h3>
                  </div>
                  <div class="box-body">
                     <img class="img-responsive pad" src="q2-summarizer-resources/histogram.png">
                     <a href="q2-summarizer-resources/histogram.png" class="btn btn-sm bg-maroon btn-flat pull-left">Open</a>
                  </div>
                </div>
              </div>
              <div class="col-md-6">
                <div class="box box-solid">
                  <div class="box-header with-border">
                    <h3 class="box-title">OTU Rank Abundance</h3>
                  </div>
                  <div class="box-body">
                    <!-- OTU Rank Abundance -->
                     <img class="img-responsive pad" src="q2-summarizer-resources/rank_abundance.png">
                  </div>
                </div>
              </div>
            </div>

            <!-- Sample depths table -->
            <div class="row">
              <div class="col-md-12">
                <div class="box box-solid">
                  <div class="box-header with-border">
                    <h3 class="box-title">Sampling Depth Histogram</h3>
                  </div>
                  <div class="box-body">
                     <p>table of depths </p>
                  </div>
                </div>
              </div>
            </div>
          </section>
        </div>
      </div>
      <footer class="main-footer">
        <div class="container">
          <div class="pull-right hidden-xs">
            <strong><a href="http://almsaeedstudio.com">Made with Admin LTE 2.3.6</a></strong>
          </div>
        </div>
      </footer>
    </div>

    <!-- jQuery 2.2.3 -->
    <script src="q2-summarizer-resources/dist/plugins/jQuery/jquery-2.2.3.min.js"></script>
    <!-- Bootstrap 3.3.6 -->
    <script src="q2-summarizer-resources/dist/bootstrap/js/bootstrap.min.js"></script>
    <!-- AdminLTE App -->
    <script src="q2-summarizer-resources/dist/js/app.min.js"></script>
    </body>
    </html>
    '''

    # Copy support files to directory
    current_file_path = abspath(__file__)
    current_dir_path = dirname(current_file_path)
    copy_tree(join(current_dir_path, "dist/"), join(output, "dist/"))

    # Save histogram file to disk
    sns.set_style("whitegrid")
    sns_plot = sns.distplot(depths, kde=False, rug=True, hist=True, color='#3498db')
    sns_plot.set(xlabel='Sampling depths', ylabel='Frequency')
    fig = sns_plot.get_figure()
    fig.savefig(join(output, "histogram.png"), dpi=300)

    # Save rank abundance to directory

    # Write HTML file to disk
    file = open(join(output_dir, 'index.html') ,'w')
    file.write(html)
    file.close()

    # Return nothing
    return None
