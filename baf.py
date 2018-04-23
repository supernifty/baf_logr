#!/usr/bin/env python
'''
  generate interactive baf and logr plots
'''

import argparse
import logging
import sys

import cyvcf2
import plotly

import plotly.offline as offline
import plotly.graph_objs as go

##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">

def plot_chrom(chrom, xs, ys, target):
  logging.info('plotting {} values for chromosome {}...'.format(len(xs), chrom))
  data = [go.Scatter(x=xs, y=ys, mode='markers', marker = {'size': 1})]
  layout = {'title': 'Chromosome {}'.format(chrom), 'yaxis': {'range': [0.0, 1.0]}}
  fig = {'data': data, 'layout': layout}
  offline.plot(fig, filename='{}.{}.html'.format(target, chrom), auto_open=False)
  logging.info('plotting {} values for chromosome {}: done'.format(len(xs), chrom))

def main(baf_target, logr_target, min_qual=None):
  logging.info('reading from stdin...')
  xs = []
  ys = []
  chrom = None
  skipped = 0
  for idx, variant in enumerate(cyvcf2.VCF('-')):
    if chrom is not None and variant.CHROM != chrom and len(xs) > 0:
      plot_chrom(chrom, xs, ys, baf_target)
      xs = []
      ys = []
    chrom = variant.CHROM
    if min_qual is not None:
      if float(variant.format('RBQ')) < min_qual or float(variant.format('ABQ')) < min_qual:
        skipped += 1
        continue
    baf = float(variant.format('AD')) / (float(variant.format('AD')) + float(variant.format('RD')))
    xs.append(variant.POS)
    ys.append(baf)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='generate baf and logr plots')
  parser.add_argument('--baf', default='baf')
  parser.add_argument('--logr', default='logr')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(args.baf, args.logr)
