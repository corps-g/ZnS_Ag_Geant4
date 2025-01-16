from nice_plot import *
from numpy.random import rand

# ZnS(Ag) grain placement in the scintillation volume
# Hornyak button

# dimension in cm
inch = 2.54
lenx = 5. / 8. * inch
leny = 7. / 64. * inch
lenz = 1. * inch
rZSA = 20.0e-4
dZSA = 2.0 * rZSA

wr = 0.05
rhoZSA = 4.09
rhoPMMA = 1.19

vol = lenx * leny * lenz
volZSA = wr * rhoPMMA * vol / ((1.0 - wr) * rhoZSA + wr * rhoPMMA)

grain_vol = 4. / 3. * np.pi * rZSA ** 3
# total number of grains
nZSA = volZSA / grain_vol

# fnx / lenx = fny / leny = fnz / lenz = k
k = (nZSA / vol) ** (1. / 3.)

# number of cells in z axis
fnz = k * lenz

""" The number of layers in the z axis is an int, nz.
In each z layer, the number of y layers is sampled to fny.
In each z, y layer, the number of x layers is sampled to fnx.
The choose of length axis (z) is because it is the longest dimension, which has 
the largest number of layers.
Two advantages:
1) allows enough layer number to sample fny and fnx
2) effect of the round of (fnz -> nz) is smallest.
   loss_n = (fnz - nz) * fnx * fny
   The lost grains are added into xy payer anyway.
"""

# number of z layers
nz = int(fnz)
diffz = fnz - nz
nz1 = nz + 1
print 'nz = %i' % nz
# length of a z layer
clz = lenz / nz
clz1 = lenz / nz1
assert clz1 > dZSA
rsz = clz - dZSA
rsz1 = clz1 - dZSA

# need to precompute the dz to determine remainz
dz = []
usedz = 0.0
for z in range(nz):
    if rand() < diffz:
        dz.append(1)
        usedz += clz1
    else:
        dz.append(0)
        usedz += clz
remainz = lenz - usedz
# divide the paitial cell into nz cells
unit_remain = remainz / nz
clz += unit_remain
clz1 += unit_remain

print 'clz, clz1, remainz = ', clz, clz1, remainz
# sampling length
sample_length_z = clz - dZSA
# number of grains in a z layer.
fxy = nZSA / nz

# Distribute the nxy into x and y proportional to the lenx, leny
# fnx / lenx = fny / leny = kxy
# kxy^2 = fnx * fny / (lenx * leny) = fxy / (lenx * leny)
# fnx = kxy * lenx
# fny = kxy * leny
kxy = (fxy / (lenx * leny)) ** 0.5
fnx = kxy * lenx
fny = kxy * leny
print 'fnx, fny = ', fnx, fny

# int number of grains in each direction
nx = int(fnx)
diffx = fnx - nx # probability of number = nx1
nx1 = nx + 1
# cell length
clx = lenx / nx
clx1 = lenx / nx1
assert clx1 > dZSA
# space of a x cell can be randomly subdivided to 
# 1) be the offset
# 2) sampling the center of a sphere in a normal cell
rsx = clx - dZSA
rsx1 = clx1 - dZSA


ny = int(fny)
ny1 = ny + 1
diffy = fny - ny
cly = leny / ny
cly1 = leny / ny1
assert cly1 > dZSA
rsy = cly - dZSA
rsy1 = cly1 - dZSA

print """clx = %.6e, clx1 = %.6e, 
cly = %.6e, cly1 = %.6e,
clz = %.6e (um)""" \
% (clx * 1e4, clx1 * 1e4, cly * 1e4, cly1 * 1e4, clz * 1e4)

edgez = -0.5 * lenz
edgey = -0.5 * leny
edgex = -0.5 * lenx

startz = edgez + rZSA
cnt = 0
for z in range(nz):
    # sample the number of y cells in this z cell
    if rand() < diffy:
        # ny + 1 cells
        numy = ny1
        # normal cell width
        ncy = cly1
        sample_length_y = rsy1
        # divide a normal cell into an offset and a subcell
        offsety = rand() * rsy1
        # sub cell length to place a grain
        scy = ncy - offsety 
    else:
        numy = ny
        # normal cell width
        ncy = cly
        sample_length_y = rsy
        # divide a cell into offset and subcell
        # a grain in the subcell
        offsety = rand() * rsy
        scy = ncy - offsety
    # place the offset at the head and at the tail repeatly
    # heady = 1: offset at the head
    # heady = 0: offset at the tail
    heady = z % 2
    starty = edgey + rZSA
    for y in range(numy + 1):
        if y == 0:
            if heady:
                # offset, do not place a grain
                starty += offsety
                continue
            else:
                # a sub cell
                dy = scy
                sly = scy - dZSA
        elif y < numy:
            # a normal cell
            sly = sample_length_y
            dy = ncy
        else:
            if heady:
                dy = scy
                sly = scy - dZSA 
            else:
                continue
                
        # sample the number of grains along x axis in this y layer
        if rand() < diffx:
            numx = nx1
            ncx = clx1
            sample_length_x = rsx1
            offsetx = rand() * rsx1
            # length of the subcell placing a grain
            scx = ncx - offsetx
        else:
            numx = nx
            ncx = clx
            sample_length_x = rsx
            offsetx = rand() * rsx
            scx = ncx - offsetx
        # place the offset at the beginning and at the end repeatly
        # headx = 1: offset at the head
        # headx = 0: offset at the tail
        headx = y % 2
        startx = edgex + rZSA
        for x in range(numx + 1):
            if x == 0:
                if headx:
                    # offset cell
                    startx += offsetx
                    continue
                else:
                    dx = scx
                    slx = scx - dZSA
            elif x < numx:
                # a normal cell
                slx = sample_length_x
                dx = ncx
            else:
                # the last cell
                if headx:
                    # the last cell is a sub cell
                    dx = scx
                    slx = scx - dZSA
                else:
                    # offset
                    continue
            cnt += 1
            # place a grain using
            # startx, slx
            # starty, sly
            # startz, sample_length_z
            # e.g.,
            cx = startx + rand() * slx
            cy = starty + rand() * sly
            cz = startz + rand() * sample_length_z
            startx += dx
        starty += dy
    startz += clz
error = abs(cnt - nZSA) / nZSA
print """
Set number of grains = %.6f
Placed number of grains = %i
Relative error = %.6e
""" % (nZSA, cnt, error)

# immersion detector
print 'Immersion detector'
# dimensions in cm
lenx = 0.1255 * 2.0
leny = 0.4445 * 2.0
lenz = 5.0
wr = 0.6
 
vol = lenx * leny * lenz
volZSA = wr * rhoPMMA * vol / ((1.0 - wr) * rhoZSA + wr * rhoPMMA)
nZSA = volZSA / (4. / 3. * 3.141592653589793 * rZSA * rZSA * rZSA)
k = (nZSA / vol) ** (1. / 3.)
fnz = k * lenz
nz = int(fnz)
clz = lenz / nz

fxy = nZSA / nz
kxy = (fxy / (lenx * leny)) ** 0.5
fnx = kxy * lenx
fny = kxy * leny

clx = lenx / fnx
print 'cell length x = ', clx
pcx = lenx - clx * int(fnx)
print 'partial cell x = ', pcx

print 'nZSA = %i' % nZSA



















