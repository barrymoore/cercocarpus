#! /usr/bin/env python
import fire
import sys

class Cercocarpus():

  def echo_me(self, file):
    
    fh = None
    if (file == '-'):
      fh = sys.stdin
    else:
      fh = open(file, 'r')
      
    for line in fh: 
      print(line)

if __name__ == '__main__':
  cerc = Cercocarpus()
  fire.Fire({'echo': cerc.echo_me})
