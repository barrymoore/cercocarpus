#! /usr/bin/env python
import fire
import sys

class Cercocarpus():

  def __init__(self):
    pass
    
  def echo(self, file='-'):
    
    fh = None
    if (file == '-'):
      fh = sys.stdin
    else:
      fh = open(file, 'r')
      
    for line in fh: 
      print(line)

if __name__ == '__main__':
  cerc = Cercocarpus()
  fire.Fire({'echo': cerc.echo})
