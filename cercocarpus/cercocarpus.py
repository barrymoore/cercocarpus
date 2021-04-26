#! /usr/bin/env python
import fire
import sys

class Cercocarpus(object):

  def __init__(self):
    pass
    
  def echo(self, file='-'):
    
    fh = None
    if (file == '-'):
      fh = sys.stdin
    else:
      fh = open(file, 'r')
      
    for line in fh:
      print(line.rstrip())

  def txt2comma(self, file='-'):
    
    fh = None
    if (file == '-'):
      fh = sys.stdin
    else:
      fh = open(file, 'r')

    rows = [x.rstrip() for x in fh.readlines()]
    print(','.join(rows))

  def txt2table(self, file='-'):
    
    import pandas as pd
    from tabulate import tabulate
    
    fh = None
    if (file == '-'):
      fh = sys.stdin
    else:
      fh = open(file, 'r')

    df = pd.read_table(fh, header=None)
    print(tabulate(df.values.tolist()))

  def transpose(self, file='-'):
    
    import pandas as pd
    
    fh = None
    if (file == '-'):
      fh = sys.stdin
    else:
      fh = open(file, 'r')

    df = pd.read_table(fh, header=None).transpose()
    for index, row in df.iterrows():
      print('\t'.join([str(x) for x in row.tolist()]))

      
if __name__ == '__main__':
  cerc = Cercocarpus()
  fire.Fire({
    'echo': cerc.echo,
    'txt2comma': cerc.txt2comma,
    'txt2table': cerc.txt2table,
    'transpose': cerc.transpose,
  })
