#!/usr/bin/env python3
'''
Created on Nov 9, 2020 by C. Melnick
'''
import sys, os
import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import json
from matplotlib import colors

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
                     
                     
    group.add_argument("--hyb", dest="is_hyb",
                     default=False,
                     action="store_true",
                     help="plot hybridisation functions")
                     
    group.add_argument("--vertex", dest="is_vertex",
                     default=False,
                     action="store_true",
                     help="plot a slice of vertex_key at omega_boson = cutoff")
                     
    group.add_argument("--key", dest="key",
                     default="",
                     help="Name vertex component to plot")
                     
    group.add_argument("--cutoff", dest="cutoff",
                     default=50, type=int,
                     help="Show firs {cutoff} frequencies of {field}")
                     
    return parser

def scatter(ax, func, label):
    n = len(func)
    x = list(range(n))
    ax.plot(x,func, label = label)
    
def image(fig, ax, func, nf, nb, ib):
    
    func = np.array(func)
    func = func.reshape(nb,nf,nf)
    nr = 40
    mid = int((nf+ib)/2)
    if 0 and np.min(func) < 0 and np.max(func) > 0:
        divnorm=colors.TwoSlopeNorm(vmin=np.min(func), vcenter=0, vmax=np.max(func))
        cs = ax.imshow(func[ib,mid-nr:mid+nr,mid-nr:mid+nr], interpolation='bilinear', cmap = "coolwarm", norm=divnorm)
    else:
        cs = ax.imshow(func[ib,:,:], interpolation='bilinear')
    fig.colorbar(cs, ax=ax)

def main():
    
    parser = Arguments(add_help=True)
    opts = parser.parse_args(sys.argv[1:])
        
    obs = None
    with open(f"{opts.name}.obs.json", "r") as f:
        data = f.read()
        obs = json.loads(data)
        
    if opts.is_vertex:
        with open(f"{opts.name}.json", "r") as f:
            data = f.read()
            data = json.loads(data)
            try:
                space = opts.space.split(" ")
                if len(space)>1:
                    space = " ".join(space[:-1])
                else:
                    space = space[0]
                nb = int(data[space]["boson cutoff"])
                nf = 2*int(data[space]["fermion cutoff"])
            except:
                space = "kernels"
                nb = int(data[space]["vertex boson cutoff"])
                nf = 2*int(data[space]["vertex fermion cutoff"])
        
    field = obs[opts.space]
    for k in opts.field.split(","):
        field=field[k]
    
    if opts.field == "expansion histogram":
        fig,ax = plt.subplots(1,1)
        scatter(ax,field, "expansion histogram")
        plt.show()
        return
        
    if opts.is_vertex:
        if opts.key == "":
            for key in field.keys():
                fig,ax = plt.subplots(2,1)
                re = image(fig, ax[0], field[key]["function"]["real"], nf, nb, opts.cutoff)
                im = image(fig, ax[1], field[key]["function"]["imag"], nf, nb, opts.cutoff)
                ax[0].set_title(key)
        elif opts.key == "charge":
            for sign in [-1,1]:
                fig,ax = plt.subplots(2,1)
                key1 = "0_1_0_1"
                key2 = "0_1_2_3"
                re = image(fig, ax[0], np.array(field[key1]["function"]["real"]) + sign*np.array(field[key2]["function"]["real"]), nf, nb, opts.cutoff)
                im = image(fig, ax[1], np.array(field[key1]["function"]["imag"]) + sign*np.array(field[key2]["function"]["imag"]), nf, nb, opts.cutoff)
                ax[0].set_title(sign)
        else:
            key = opts.key
            fig,ax = plt.subplots(2,1)
            re = image(fig, ax[0], field[key]["function"]["real"], nf, nb, opts.cutoff)
            im = image(fig, ax[1], field[key]["function"]["imag"], nf, nb, opts.cutoff)
            ax[0].set_title(key)
            
        
        plt.show()
        return
    
        
    if opts.key == "":
        keys = field.keys()
    else:
        keys = [opts.key]
    
    fig,ax = plt.subplots(2,1)
    for key in keys:
        
        im = field[key]["function"]["imag"]
        re = field[key]["function"]["real"]
        
        cutoff = min(opts.cutoff, len(im))
        
        scatter(ax[0],im[:cutoff], key)
        ax[0].set_ylabel(f"Im[{opts.field}]")
        
        scatter(ax[1],re[:cutoff], key)
        ax[1].set_ylabel(f"Re[{opts.field}]")
        ax[1].set_xlabel(f"n")
        ax[0].set_title(key)

    plt.legend()
    plt.show()
    
    
    
if __name__ == '__main__':
    main()
    
