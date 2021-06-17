#!/usr/bin/env python3
'''
Created on Nov 9, 2020 by C. Melnick
'''
import sys, os
import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import json

def Arguments(add_help=False):
    
    parser = ArgumentParser('ComCTQMC', add_help=add_help,
                            formatter_class=ArgumentDefaultsHelpFormatter)    
    
    group = parser.add_argument_group("Plot controls")

    group.add_argument("--space", dest="space",
                     default="partition",
                     help="configuration space from which to plot an observable")
                    
    group.add_argument("--field", dest="field",
                     default="self-energy",
                     help="Observable to plot (from space)")
                     
    group.add_argument("--run-name", dest="name",
                     default="params",
                     help="Name of run. (input file will be named [this].params)")
                     
                     
    group.add_argument("--cutoff", dest="cutoff",
                     default=50, type=int,
                     help="Show firs {cutoff} frequencies of {field}")
                     
    return parser

def scatter(ax, func):
    
    n = len(func)
    x = list(range(n))
    ax.scatter(x,func)

def main():
    
    parser = Arguments(add_help=True)
    opts = parser.parse_args(sys.argv[1:])
    
    obs = None
    with open(f"{opts.name}.obs.json", "r") as f:
        data = f.read()
        obs = json.loads(data)
        
    field = obs[opts.space][opts.field]
    
    if opts.field == "expansion histogram":
        fig,ax = plt.subplots(1,1)
        scatter(ax,field)
        plt.show()
        return
    
    for key in field.keys():
        fig,ax = plt.subplots(2,1)
        im = field[key]["function"]["imag"]
        re = field[key]["function"]["real"]
        
        cutoff = min(opts.cutoff, len(im))
        
        scatter(ax[0],im[:cutoff])
        ax[0].set_ylabel(f"Im[{opts.field}]")
        
        scatter(ax[1],re[:cutoff])
        ax[1].set_ylabel(f"Re[{opts.field}]")
        ax[1].set_xlabel(f"n")
        
    plt.show()
    
    
    
if __name__ == '__main__':
    main()
    
