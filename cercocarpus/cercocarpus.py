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

      
if __name__ == '__main__':
  cerc = Cercocarpus()
  fire.Fire({'echo': cerc.echo,
             'txt2comma': cerc.txt2comma})
