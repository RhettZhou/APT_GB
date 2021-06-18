#-----------------------------------------------------------
# File mesh.py
#-----------------------------------------------------------
# create mesh function wrappers

import numpy as np
import bmesh
import os
from math import floor, sqrt
from collections import Counter
from bpy.props import *
from mathutils import Vector, Matrix
from scipy.spatial import Delaunay, ConvexHull

# A function to cut the wedge with a specific mesh width
def cut_edge(me, bm, select_edge, mesh_width, epsilon):
    cut_num = floor(select_edge/mesh_width)
    for edge in bm.edges:
        edge.select = False
        if abs(round(edge.calc_length(),4)-select_edge)<=epsilon:
            edge.select = True
    selected_edges = [edge for edge in bm.edges if edge.select]
    bmesh.ops.subdivide_edges(bm, edges =selected_edges, cuts = cut_num)
    bmesh.update_edit_mesh(me)

def cal_length(v1,v2):
    v = []
    v.append([v2[0] - v1[0],v2[1] - v1[1],v2[2] - v1[2]])
    v_length = sqrt(v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2])
    return v_length

def cal_unit(v):
    v_length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    v = v/v_length
    return v

def cal_nor(v1,v2):
    nor = np.cross(v1,v2)
    nor = cal_unit(nor)
    return nor

def cal_hull(verts,incre):
    vect1 = cal_unit(verts[1] - verts[0])
    vect2 = cal_unit(verts[2] - verts[1])
    gb_h = []                                                    # expand in plane
    gb_h.append([verts[0][0] - incre*vect1[0] - incre*vect2[0], verts[0][1] - incre*vect1[1] - incre*vect2[1], verts[0][2] - incre*vect1[2] - incre*vect2[2]])
    gb_h.append([verts[1][0] + incre*vect1[0] - incre*vect2[0], verts[1][1] + incre*vect1[1] - incre*vect2[1], verts[1][2] + incre*vect1[2] - incre*vect2[2]])
    gb_h.append([verts[2][0] + incre*vect1[0] + incre*vect2[0], verts[2][1] + incre*vect1[1] + incre*vect2[1], verts[2][2] + incre*vect1[2] + incre*vect2[2]])
    gb_h.append([verts[3][0] - incre*vect1[0] + incre*vect2[0], verts[3][1] - incre*vect1[1] + incre*vect2[1], verts[3][2] - incre*vect1[2] + incre*vect2[2]])
    nor = cal_nor(vect1, vect2)
    gb_verts = []                                                # expand out of plane
    for i in range(4):
        gb_verts.append([gb_h[i][0] - incre*nor[0], gb_h[i][1] - incre*nor[1], gb_h[i][2] - incre*nor[2]])
        gb_verts.append([gb_h[i][0] + incre*nor[0], gb_h[i][1] + incre*nor[1], gb_h[i][2] + incre*nor[2]])
    gb_hull = Delaunay(gb_verts)
    return gb_hull  

def in_hull(p, hull):
    return hull.find_simplex(p)>=0

def obj_rename(name):
    filename = os.path.splitext(os.path.basename(name))
    items = filename[0].split('__')
    item0 = items[0].split('_')[0]
    item1 = items[1].split('_')[0]
    filename_short = [item0 + '_' + item1]
    return filename_short

def obj_filepath(name):
    filepath = os.path.dirname(name)
    return filepath

def obj_absname(importName,export_obj_name,quad,target,radius):
    folder = obj_filepath(importName)
    obj = export_obj_name

    export_obj_absname = [folder + '\\' + obj + '_QLen' + str(quad) + '_' + str(target) + '_R' + str(radius) + '.obj']

    return export_obj_absname[0] 




