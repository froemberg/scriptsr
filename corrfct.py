#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2008-2012 Felix Höfling
#

import sys, math
from numpy import *
import argparse
import os
from os import makedirs, path


desc = """
Usage: correlation.py [options] [data_file]

Calculate autocorrelation of a time series,
read <stdin> if no data file specified
"""
parser = argparse.ArgumentParser(description=desc)

parser.add_argument('input', metavar='INPUT', nargs='?', default='/dev/stdin', help='input file')

parser.add_argument('-o', dest='output_file', default='/dev/stdout',
                    help='write result to FILE', metavar='FILE')

parser.add_argument('-c', '--column', type=int, nargs='+', default=(0,), dest='col', metavar='<n>',
                   help='specify column(s) of data file to be used')

parser.add_argument('--acf', action='store_const', const='acf', dest='action', default='acf',
                    help='calculate autocorrelation function <x(t)x(0)> / <x^2>')

parser.add_argument('--msd', action='store_const', const='msd', dest='action',
                    help='calculate mean-square displacement <[x(t)-x(0)]^2>')

parser.add_argument('--psd', action='store_const', const='psd', dest='action',
                    help='calculate power spectral density <|x(omega)|^2> / T')

parser.add_argument('--increments', action='store_true', dest='increments', default=False,
                    help='use increments between data points')

parser.add_argument('--fft', action='store_true', dest='fft', default=False,
                    help='use FFT-based algorithm. Caveat: implies periodic continuation of data set')

parser.add_argument('--resolution', type=float, default=1, metavar='<x>',
                    help='time resolution of sampling')

parser.add_argument('--block-size', type=int, dest='block_size', default=1, metavar='<n>',
                    help='block size for coarse-graining output on logarithmic time axis')

parser.add_argument('--block-factor', type=int, dest='block_factor', default=1, metavar='<x>',
                    help='coarse-graining factor')

parser.add_argument('--input-size', type=int, dest='input_size', metavar='<n>',
                    help='delimit size of input data, use more efficient algorithm')

parser.add_argument('--revreaddy-input', action='store_true', dest='revreaddy_input', default=False,
                    help='read HDF5 input file from revReaDDy 0.2 software')

parser.add_argument('-a', '--append', action='store_true', dest='append', default=False,
                    help='append to output file')

parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', default=False,
                    help='print verbose messages')

args = parser.parse_args()










data_dir = os.path.join(path.dirname(os.path.realpath(__file__)), "..", "data/sim05q")

output_path = os.path.join(data_dir, "corrdata")
if not os.path.exists(output_path):
    makedirs(output_path)

infile_path = os.path.join(data_dir, args.input)


# read data from input file
if args.verbose:
    if len(args.col) == 1:
        colstring = ' {0}'.format(args.col[0])
    else:
        colstring = 's ' + ', '.join(['{0}'.format(x) for x in args.col])
    print 'Reading data from column{0} of file {1:s}'.format(colstring, infile_path)

if not args.revreaddy_input:
    if not args.input_size:
        # read complete text file using NumPy's default function
        data = loadtxt(infile_path, usecols=args.col, dtype=float32)
    else:
        # read given number of lines from a text file using Python's I/O functions,
        # this is expected to be faster than loadtxt()
        from itertools import ifilterfalse
        def ignoreline(s):
            s = s.strip()
            return s.startswith('#') or len(s) == 0

        data = empty((args.input_size, len(args.col)), dtype=float32)
        col = array(args.col)
        with open(file_path) as f:
            for i,line in enumerate(ifilterfalse(ignoreline, f)):
                data[i]=array(line.split(), dtype=data.dtype)[col]
                if i + 1 >= data.shape[0]:
                    break
else:
    # read HDF5 file produced by the revReaDDy 0.2 software
    import h5py

    with h5py.File(infile_path, 'r') as f:
        positions = f['positions']
        size = positions.shape[0]
        if args.input_size:
            size = min(args.input_size, size)

        if args.verbose:
            print 'Reading {0} data points from HDF5 file'.format(size)

        data = positions[:size, 0, args.col[0]]   # particle #0, args.col is the Cartesian component

if args.increments:
    data = data[1:] - data[:-1]
    if args.verbose:
        'using increments between data points'

N = data.shape[0]



# calculate correlations




if args.verbose:
    print 'data set contains %d values' % N
    print 'mean={0}, std deviation={1}'.format(mean(data, axis=0, dtype=float64), std(data, axis=0, dtype=float64))

    print 'Calculating', \
          (args.action=='msd' and 'mean-square displacement...') \
            or (args.action=='acf' and 'normalized autocorrelation...') \
            or (args.action=='psd' and 'power spectral density...')

if args.fft or args.action=='psd':
    ft_data = fft.rfft(data, axis=0)
    spectrum = (ft_data * ft_data.conj()) / N

if args.action=='msd':
    if args.fft:
        corr = fft.irfft(spectrum, axis=0)[:N/2]
        msd = 2*(corr[0]-corr)
    else:
        msd = empty((N-1,) + tuple(data.shape[1:]))
        for n in range(N-1):
            A = data[n:]
            B = data[:N-n]
            tmp = A - B
            msd[n] = average(tmp*tmp, axis=0)
    result = msd

elif args.action=='acf':
    if args.fft:
        corr = fft.irfft(spectrum, axis=0)[:N/2]
    else:
        corr = empty(N-1)
        for n in range(N-1):
            A = data[n:]
            B = data[:N-n]
            corr[n] = average(A*B, axis=0)
    result = corr

elif args.action=="psd":
    spectrum = spectrum.real
    result = spectrum * args.resolution

# write results to output file
file = open(os.path.join(output_path, args.output_file + ".dat"), (args.append and 'a') or 'w')



if args.verbose:
    print '{mode:s} results to file {file:s}'.format(
        mode=(args.append and "Appending") or "Writing",
        file=args.output_file
    )

print >>file, (args.action=='msd' and '# time MSD') \
  or (args.action=='acf' and '# time ACF') \
  or (args.action=='psd' and '# frequency PSD')

# coarse-grain result on logarithmic time axis
a = 0
n = 0
m = 1
while n+m <= result.shape[0]:
    avg = sum(average(result[n:n+m], axis=0))
    err = sqrt(sum(var(result[n:n+m], axis=0)))
    # err = m>1 and std(result[n:n+m], axis=0)/sqrt(m-1) or 0
    if not args.action == 'psd':
        x = args.resolution * (n+.5*(m-1))  # time grid
    else:
        x = (n+.5*(m-1)) * 2*pi / (N * args.resolution)  # frequency grid (ω)

    print >>file, '{0:g}\t{1:g}\t{2:g}'.format(x, avg, err)
    n += m
    a += 1
    if a % args.block_size == 0:
        m *= args.block_factor

print >>file, '\n'

