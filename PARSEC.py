from __future__ import division
from math import *
import numpy as np
from ofblockmeshdicthelper import *
import fileinput
import sys
import os
import matplotlib.pyplot as plt


# ========================================================================================
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def pcoef(
        xte, yte, rle,
        x_cre, y_cre, d2ydx2_cre, th_cre,
        surface):
    """evaluate the PARSEC coefficients"""
    # Initialize coefficients
    coef = np.zeros(6)

    # 1st coefficient depends on surface (pressure or suction)
    if surface.startswith('p'):
        coef[0] = -sqrt(2 * rle)
    else:
        coef[0] = sqrt(2 * rle)

    # Form system of equations
    A = np.array([
        [xte ** 1.5, xte ** 2.5, xte ** 3.5, xte ** 4.5, xte ** 5.5],
        [x_cre ** 1.5, x_cre ** 2.5, x_cre ** 3.5, x_cre ** 4.5,
         x_cre ** 5.5],
        [1.5 * sqrt(xte), 2.5 * xte ** 1.5, 3.5 * xte ** 2.5,
         4.5 * xte ** 3.5, 5.5 * xte ** 4.5],
        [1.5 * sqrt(x_cre), 2.5 * x_cre ** 1.5, 3.5 * x_cre ** 2.5,
         4.5 * x_cre ** 3.5, 5.5 * x_cre ** 4.5],
        [0.75 * (1 / sqrt(x_cre)), 3.75 * sqrt(x_cre), 8.75 * x_cre ** 1.5,
         15.75 * x_cre ** 2.5, 24.75 * x_cre ** 3.5]
    ])

    B = np.array([
        [yte - coef[0] * sqrt(xte)],
        [y_cre - coef[0] * sqrt(x_cre)],
        [tan(th_cre * pi / 180) - 0.5 * coef[0] * (1 / sqrt(xte))],
        [-0.5 * coef[0] * (1 / sqrt(x_cre))],
        [d2ydx2_cre + 0.25 * coef[0] * x_cre ** (-1.5)]
    ])

    # Solve system of linear equations
    X = np.linalg.solve(A, B)

    # Gather all coefficients
    coef[1:6] = X[0:5, 0]

    # Return coefficients
    return coef


# ========================================================================================
# Sample coefficients  NACA0012
'''
leading edge radius (r_le)
pressure and suction surface crest locations (x_pre, y_pre, x_suc, y_suc)
curvatures at the pressure and suction surface crest locations (d2y/dx2_pre, d2y/dx2_suc)
trailing edge coordinates (y_te)
trailing edge angles between the pressure and suction surface and the horizontal axis (th_pre, th_suc)
'''
rle_suc = .014927
rle_pre = .014181
x_suc = .29866
y_suc = .059404
x_pre = .29962
y_pre = -0.059632
d2ydx2_suc = -.42399
d2ydx2_pre = 0.445281
yte = 0
th_suc = -7.672047
th_pre = 7.59506

xte = 1.0  # Trailing edge x position
# Evaluate pressure (lower) surface coefficients
cf_pre = pcoef(xte, yte, rle_pre,
               x_pre, y_pre, d2ydx2_pre, th_pre,
               'pre')
# Evaluate suction (upper) surface coefficients
cf_suc = pcoef(xte, yte, rle_suc,
               x_suc, y_suc, d2ydx2_suc, th_suc,
               'suc')
x = (1 - np.cos(np.linspace(0, 1, np.ceil(1e3)) * np.pi)) / 2
uppery = np.array([0] * len(x))
lowery = np.array([0] * len(x))
for i in range(1, 7):
    uppery = uppery + cf_suc[i - 1] * x ** (i - 1 / 2)
    lowery = lowery + cf_pre[i - 1] * x ** (i - 1 / 2)

# ========================================================================================
# Plot this airfoil
# plt.plot(x, uppery, 'o--', markerfacecolor='red')
# plt.plot(x, lowery, 'o--', markerfacecolor='blue')
# plt.gca().axis('equal')
# plt.show()
# ========================================================================================
# create blockmesh
# prepare ofblockmeshdicthelper.BlockMeshDict instance to
# gather vertices, blocks, faces and boundaries.
bmd = BlockMeshDict()
length_z = xte / 100
L1 = 15
L2 = 15
number = round(x_suc / 1 * len(x))
x10 = x[number]
y10 = uppery[number]
x11 = x[number]
y11 = lowery[number]


ymax=max(uppery)
ymin=min(lowery)
A=(ymax-ymin)*length_z
bmd.set_metric('m')
basevs = [Vertex(x10 - L1, 0, 0, 'v0'),
          Vertex(x11, -L1, 0, 'v1'),
          Vertex(1, -L1, 0, 'v2'),
          Vertex(1 + L2, -L1, 0, 'v3'),
          Vertex(1 + L2, 0, 0, 'v4'),
          Vertex(1 + L2, L1, 0, 'v5'),
          Vertex(1, L1, 0, 'v6'),
          Vertex(x10, L1, 0, 'v7'),
          Vertex(0, 0, 0, 'v8'),
          Vertex(1, 0, 0, 'v9'),
          Vertex(x10, y10, 0, 'v10'),
          Vertex(x11, y11, 0, 'v11'),
          ]

for v in basevs:
    bmd.add_vertex(v.x, v.y, v.z, v.name + '-z')
    bmd.add_vertex(v.x, v.y, v.z + length_z, v.name + '+z')

z = length_z / 2
for i in ['+z', '-z']:
    bmd.add_arcedge(('v0' + i, 'v7' + i), 'arc-h' + i,
                    Vertex(x10 - L1 * sin(pi / 6), L1 * cos(pi / 6), length_z / 2 + eval(i), 'arch' + i))
    bmd.add_arcedge(('v0' + i, 'v1' + i), 'arc-l' + i,
                    Vertex(x10 - L1 * sin(pi / 6), -L1 * cos(pi / 6), length_z / 2 + eval(i), 'arcl' + i))

    upperpoints1 = []
    for j in range(1, number):
        upperpoints1.append(Point(x[j], uppery[j], length_z / 2 + eval(i)))
    lowerpoints1 = []
    for j in range(1, number):
        lowerpoints1.append(Point(x[j], lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v8' + i, 'v10' + i), 'uppercurve-1' + i, upperpoints1)
    bmd.add_splineedge(('v8' + i, 'v11' + i), 'lowercurve-1' + i, lowerpoints1)

    upperpoints2 = []
    for j in range(number + 1, len(x)):
        upperpoints2.append(Point(x[j], uppery[j], length_z / 2 + eval(i)))
    lowerpoints2 = []
    for j in range(number + 1, len(x)):
        lowerpoints2.append(Point(x[j], lowery[j], length_z / 2 + eval(i)))
    bmd.add_splineedge(('v10' + i, 'v9' + i), 'uppercurve-2' + i, upperpoints2)
    bmd.add_splineedge(('v11' + i, 'v9' + i), 'lowercurve-2' + i, lowerpoints2)


def vnamegen(x0z0, x1z0, x1z1, x0z1):
    return (x0z0 + '-z', x1z0 + '-z', x1z1 + '-z', x0z1 + '-z',
            x0z0 + '+z', x1z0 + '+z', x1z1 + '+z', x0z1 + '+z')


# create blocks
# airfoil to far field
yCells =   62*3
yGrading = 200
# https://turbmodels.larc.nasa.gov/naca0012numerics_grids.html
xCell_total = int(yCells*2/(1.49))

# downstream
xDCells = int(xCell_total*0.3)
xDGrading = 100

xCell_airfoil = int(xCell_total*0.3)
# upstream
xUCells = int(3*x10 * xCell_airfoil)
xUleading = 600
#xUleading = 400
leadGrading1= 0.5
leadGrading2= 0.6
# middle
xMCells = xCell_airfoil - xUCells
b0 = bmd.add_hexblock(vnamegen('v0', 'v1', 'v11', 'v8'), (xUCells, yCells, 1), 'b0',
                      EdgeGrading(1 / leadGrading1, 1 / leadGrading1, 1 / leadGrading1, 1 / leadGrading1, 1 / xUleading,
                                  1 / yGrading, 1 / yGrading,
                                  1 / xUleading, 1, 1, 1, 1))
b1 = bmd.add_hexblock(vnamegen('v1', 'v2', 'v9', 'v11'), (xMCells, yCells, 1), 'b1', SimpleGrading(leadGrading2, 1 / yGrading, 1))
b2 = bmd.add_hexblock(vnamegen('v2', 'v3', 'v4', 'v9'), (xDCells, yCells, 1), 'b2',
                      SimpleGrading(xDGrading, 1 / yGrading, 1))
b3 = bmd.add_hexblock(vnamegen('v9', 'v4', 'v5', 'v6'), (xDCells, yCells, 1), 'b3',
                      SimpleGrading(xDGrading, yGrading, 1))
b4 = bmd.add_hexblock(vnamegen('v10', 'v9', 'v6', 'v7'), (xMCells, yCells, 1), 'b4', SimpleGrading(leadGrading2, yGrading, 1))
b5 = bmd.add_hexblock(vnamegen('v0', 'v8', 'v10', 'v7'), (yCells, xUCells, 1), 'b5',
                      EdgeGrading(1 / xUleading, 1 / yGrading, 1 / yGrading, 1 / xUleading, 1 / leadGrading1,
                                  1 / leadGrading1, 1 / leadGrading1,
                                  1 / leadGrading1, 1, 1, 1, 1))

# ========================================================================================
# face element of block can be generated by Block.face method
bmd.add_boundary('patch', 'inlet-l', [b0.face('ym')])
bmd.add_boundary('wall', 'airfoil-l1', [b0.face('yp')])
bmd.add_boundary('empty', 'empty-1', [b0.face('zm')])
bmd.add_boundary('empty', 'empty-2', [b0.face('zp')])

bmd.add_boundary('wall', 'wall-l1', [b1.face('ym')])
bmd.add_boundary('wall', 'airfoil-l2', [b1.face('yp')])
bmd.add_boundary('empty', 'empty-3', [b1.face('zm')])
bmd.add_boundary('empty', 'empty-4', [b1.face('zp')])

bmd.add_boundary('wall', 'wall-l2', [b2.face('ym')])
bmd.add_boundary('patch', 'outlet-l', [b2.face('xp')])
bmd.add_boundary('empty', 'empty-5', [b2.face('zm')])
bmd.add_boundary('empty', 'empty-6', [b2.face('zp')])

bmd.add_boundary('patch', 'outlet-h', [b3.face('xp')])
bmd.add_boundary('wall', 'wall-h2', [b3.face('yp')])
bmd.add_boundary('empty', 'empty-7', [b3.face('zm')])
bmd.add_boundary('empty', 'empty-8', [b3.face('zp')])

bmd.add_boundary('wall', 'airfoil-h2', [b4.face('ym')])
bmd.add_boundary('wall', 'wall-h1', [b4.face('yp')])
bmd.add_boundary('empty', 'empty-9', [b4.face('zm')])
bmd.add_boundary('empty', 'empty-10', [b4.face('zp')])

bmd.add_boundary('patch', 'inlet-h', [b5.face('xm')])
bmd.add_boundary('wall', 'airfoil-h1', [b5.face('xp')])
bmd.add_boundary('empty', 'empty-11', [b5.face('zm')])
bmd.add_boundary('empty', 'empty-12', [b5.face('zp')])

# ========================================================================================
# prepare for output
bmd.assign_vertexid()
# output
f = open('blockMeshDict', 'wb')
f.write(bmd.format().encode())
f.close()
