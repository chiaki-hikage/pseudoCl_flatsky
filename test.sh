#!/bin/sh

CDIR=`pwd`
OUTPUT=$CDIR

### number of cores for parallel computing
ppn=4

### Noise Add ?  [y/n]

addnoi=n

### periodic boundary condition [y/n]

pcond=n

### intial seed to generate random number

irun=0

###  calc mode couping matrix due to a square boundary
###  (no need if periodic boundary condition is assumed)

if [ $pcond = "y" ]; then
mat_square='none'
else
mat_square=$OUTPUT'/modemat_square_para.dat'
mpirun -np $ppn ./modemat_square_para $mat_square
fi

###  calculation start

sumpow=$OUTPUT'/sumpow.dat'

### make Gaussian shear data

sfile=$OUTPUT'/initshear.dat'
ipower=$OUTPUT'/inputpow.dat'
./sheardist $sfile $ipower $irun $pcond

###  make simulated mask from the shear data

mfile=$OUTPUT'/mask.dat'
ffile=$OUTPUT'/maskflag.dat'
./starmask $sfile $mfile $ffile

###  calculate mode-coupling matrix due to the simulated mask

mat_inside=$OUTPUT'/modemat_inside_para.dat'
mpirun -np $ppn ./modemat_inside_para $mfile $mfile $mat_inside

###  calculate convolved shear power spectrum

cpower=$OUTPUT'/convpow.dat'
./deconv $ipower 'none' $cpower $mat_square $mat_inside 'c'

###  calculate masked shear power spectrum

mpower=$OUTPUT'/maskpow.dat'
./shearpow $sfile $mpower $addnoi $ffile 'n' $irun

###  calculate noise power spectrum by rotating shear at random

if [ $addnoi = "y" ]; then
npower=$OUTPUT'/noisepow.dat'
./shearpow $sfile $npower $addnoi $ffile 'y' $irun
else
npower='none'
fi

###  deconvolve masked shear power spectrum

dpower=$OUTPUT'/deconvpow.dat'
./deconv $mpower $npower $dpower $mat_square $mat_inside 'd'

###  gather powers into $sumpow

./collect $ipower $dpower $cpower $mpower $sumpow
