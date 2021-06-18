# =============================================================================
# (C) Copyright 2014
# Australian Centre for Microscopy & Microanalysis
# The University of Sydney
# =============================================================================
# File:   blend/space.py
# Date:   2014-07-01
# Author: Varvara Efremova
#         Modified by Rhett Xuyang Zhou 20210521
# Description:
# Blender API wrapper for scene/world/etc general operations
# =============================================================================

import bpy
import numpy as np
import bmesh
from math import floor, sqrt
from collections import Counter
from bpy.props import *
from mathutils import Vector

from scipy.spatial import Delaunay, ConvexHull


from .aptread import APTloader
from . import blend, mesh

# === Operator execute functions ===
# === Operator execute functions ===
def scale_child(self, context):
    obj = context.object
    child = obj.children[0]

    # select child and deselect obj
    blend.object.select(obj, False)
    blend.object.select(child, True)
    return blend.object.selected_resize()

def position_active_camera_on(self, context):
    blend.space.camera_position_on()
    return {'FINISHED'}

def position_active_camera_off(self, context):
    blend.space.camera_position_off()
    return {'FINISHED'}

def make_camera_active(self, context):
    cam = context.object
    blend.space.camera_set_active(cam)
    return {'FINISHED'}

def add_bounding_box(self, context):
    """Calculate and add a bounding box to the current data"""
    props = context.scene.pos_panel_props
    padding = self.padding

    # FIXME don't load this again!!! save as global var for now?
    data = APTloader.ReadAPTData(props.pos_filename, props.rng_filename)

    pointlist = data.xyz
    xyzmax = np.amax(pointlist, axis=0) # max locations in data
    xyzmin = np.amin(pointlist, axis=0) # min locations in data

    # add padding to extremal xyz coords
    xyzmax += padding
    xyzmin -= padding

    # define vertices of cube, clockwise, starting from bottom
    coords = [[-1, -1, -1],
             [-1,  1, -1],
             [ 1,  1, -1],
             [ 1, -1, -1],
             [-1, -1,  1],
             [-1,  1,  1],
             [ 1,  1,  1],
             [ 1, -1,  1]]
    coords = np.array(coords)
    # numpy where is awesome
    # (replaces -ves with respective min values, +ves with max values)
    verts = np.where(coords > 0, xyzmax, xyzmin)

    # create bounding box and draw as wireframe
    name = "Bound"
    origin = (0.0, 0.0, 0.0)
    bound = blend.object.cube_add_from_verts(name, origin, verts)
    boundmat = blend.material.surface_add(name+"_mat", shadeless=True)
    blend.material.set(bound, boundmat)
    blend.object.modifier_add_wireframe(bound)

    grp = blend.space.group_get("Bounds")
    if not grp:
        grp = blend.space.group_add("Bounds")
    blend.space.group_add_object(grp, bound)

    bound.datatype = 'BOUND'
    return {'FINISHED'}

def add_lamp_view(self, context):
    props = context.scene.pos_panel_props

    # FIXME temp hack
    cam, lamp = blend.space.camlamp_add_to_view(ltype='SUN')
    blend.object.delete(cam)
    context.scene.objects.active = lamp
    # end temp hack

    # proper code ...
    #drawing.add_lamp_to_view(name="Lamp", ltype='SUN', color=color)

    # Add to lamp group
    grp = blend.space.group_get("Lamps")
    if not grp:
        grp = blend.space.group_add("Lamps")
    blend.space.group_add_object(grp, lamp)
    return {'FINISHED'}

def pointcloud_add(self, context):
    props = context.scene.pos_panel_props
    color = props.ptcld_color
    emit = props.ptcld_emit

    obj = context.object

    # Get list of verts from active object
    verts = blend.object.vertices_get(obj)

    # Create pointcloud from these verts and set its parent to active obj
    ptcld = blend.object.pointcloud_add(verts, obj.name+"_pointcloud")
    blend.object.parent_set(obj, ptcld)

    # Add wireframe material
    # TODO FIXME
    mat = drawing.create_material_wire(ptcld.name+"_mat", color=color, emit=emit)
    blend.material.set(ptcld, mat)
    return {'FINISHED'}

def add_camera_view(self, context):
    cam = blend.space.camera_add_to_view(name="Camera", clip_start=0.1, clip_end=1000.0)

    # Add to camera group
    grp = blend.space.group_get("Cameras")
    if not grp:
        grp = blend.space.group_add("Cameras")
    blend.space.group_add_object(grp, cam)
    return {'FINISHED'}

def add_halo_material(self, context):
    obj = context.object

    mat = blend.material.halo_add(obj.name+"_halo", use_tex=True)
    tex = blend.material.texture_add_img(obj.name+"_tex", path=self.halo_img_path)
    blend.material.texture_add(mat, tex)
    blend.material.set(obj, mat)

    obj.vistype = 'HALO'
    return{'FINISHED'}

def remove_halo_material(self, context):
    obj = context.object
    blend.object.active_material_delete(obj)

    obj.vistype = 'NONE'
    return{'FINISHED'}

def dupli_vert(self, context):
    # Applies vertex duplication with icospheres to currently selected objects
    # select object to add mesh child to (props)
    # create desired mesh (props)

    # active object is parent
    obj = context.object

    # create child icosphere
    vertobj = blend.object.icosphere_add(name=obj.name+"_vert")

    # create material for child icosphere, link to child
    vertmat = blend.material.surface_add(name=obj.name+"_mat")
    blend.material.set(vertobj, vertmat)

    blend.object.parent_set(obj, vertobj)
    blend.object.dupli_set(obj, 'VERTS')

    obj.vistype = 'DUPLI'
    return {'FINISHED'}

def remove_duplivert(self, context):
    obj = context.object
    obj.vistype = 'NONE'
    blend.object.delete_children(obj)
    blend.object.dupli_set(obj, 'NONE')

    
    return{'FINISHED'}

def clear(self, context):
    # Clear all objects and meshes in scene
    blend.space.delete_all()
    return {'FINISHED'}

def bake(self, context):
    # Draws currently selected data
    # called by: draw_button operator
    props = context.scene.pos_panel_props
    plot_type = props.plot_type

    # FIXME temp don't do this, should be able to store ReadAPTData object in blender props??
    data = APTloader.ReadAPTData(props.pos_filename, props.rng_filename)

    if plot_type == 'ISO':
        groupname = "Isotopes"
        listfunc = "rnglist"
        getfunc = "getrng"
    elif plot_type == 'EA':
        groupname = "Elements"
        listfunc = "atomlist"
        getfunc = "getatom"
    elif plot_type == 'ION':
        groupname = "Ions"
        listfunc = "ionlist"
        getfunc = "getion"

    # populate item names and point locations
    itemlist = getattr(data, listfunc)

    namelist = []
    pointlist = []
    for ind, item in enumerate(itemlist):
        # convert item to string name
        namelist.append(str(item))
        # get points (vertices) for current item and append
        pointlist.append(getattr(data, getfunc)(ind))

    # Create group for meshes of same type
    grp = blend.space.group_add(groupname)

    # Draw all meshes in pointlist and link to group
    for name, verts in zip(namelist, pointlist):
        obj = blend.object.object_add_from_verts(verts, name, trunc=None)
        obj.datatype = 'DATA'
        blend.space.group_add_object(grp, obj)

    # centre view on created group
    blend.space.view_selected_group(groupname)
    return {'FINISHED'}

# def load_posrng(self, context):
#     # Load APT pos/rng data
#     # called by: load_posrng_button operator
#     # populates props.atomlist collection property with atoms in rng file
#     props = context.scene.pos_panel_props

#     try:
#         data = APTloader.ReadAPTData(props.pos_filename, props.rng_filename)
#         print("Loaded rng data: ", data.atomlist) # DEBUG
#         self.report({'INFO'}, "Loaded %s as POS, %s as RNG" % \
#                 (props.pos_filename, props.rng_filename))
#     except APTloader.APTReadError:
#         self.report({'ERROR'}, "Error reading pos or rng file. Double check file names.")
#         return {'CANCELLED'}

#     # separate ion names and index refs for data.getion
#     print("--- ATOMLIST", data.atomlist)
#     print("--- RNGLIST", data.rnglist)
#     print("--- IONLIST", data.ionlist)
#     return {'FINISHED'}





## ++++ Handling the obj and ply files ++++
def load_obj(self, context):
    # Load obj data
    # called by: import_obj_button operator
    props = context.scene.obj_panel_props
    
    try:
        bpy.ops.import_scene.obj(filepath=props.import_obj_filename,axis_forward='Y',axis_up='Z')
        import_obj = context.selected_objects
        import_obj_name = import_obj[0].name
        import_obj_name_short =  mesh.mesh.obj_rename(import_obj_name)
        import_obj[0].name = import_obj_name_short[0]
    except:
        return {'CANCELLED'}

    return {'FINISHED'}

def save_obj(self, context):
    # Save obj data
    # called by: export_obj_button operator
    propsO = context.scene.obj_panel_props
    propsM = context.scene.gb_mesh_panel_props

    importName = propsO.import_obj_filename
    export_obj = context.selected_objects
    export_obj_name = export_obj[0].name

    quad =propsM.quad_mesh_edge_length
    target = propsM.target_object
    radius = propsM.dcom_radius
    export_obj_absname = mesh.mesh.obj_absname(importName,export_obj_name,quad,target,radius)
    
    try:
        bpy.ops.export_scene.obj(filepath=export_obj_absname,axis_forward='Y',axis_up='Z',use_selection=True,use_mesh_modifiers=False,\
                                 use_normals=False, use_uvs=False, use_materials=False, use_vertex_groups=True)
    except:
        return {'CANCELLED'}

    return {'FINISHED'}

def load_ply(self, context):
    # Load ply data
    # called by: import_ply_button operator
    props = context.scene.obj_panel_props
    
    try:
        bpy.ops.import_mesh.ply(filepath=props.import_ply_filename)
        obj = context.object

        mat = bpy.data.materials.new(name = "Vert")
        #mat = bpy.data.materials.new()
        mat.use_nodes = True #Make so it has a node tree

        #Add the vertex color node
        vc = mat.node_tree.nodes.new('ShaderNodeVertexColor')
        #Assign its layer
        vc.layer_name = "Col"
        #Get the shader
        bsdf = mat.node_tree.nodes["Principled BSDF"]
        #Link the vertex color to the shader
        mat.node_tree.links.new( vc.outputs[0], bsdf.inputs[0] )
        #bpy.context.object.active_material.use_shadeless = True

        blend.material.set(obj, mat)
    except:
        return {'CANCELLED'}

    return {'FINISHED'}




## ++++ The following is for GBMesh ++++
def quad_mesh(self, context):
    print("")
    print("Quad mesh execution")
    
    mesh_width = bpy.context.scene.gb_mesh_panel_props.quad_mesh_edge_length            # Set parameter, mesh width
    ini_all_edges = []                                              # Define a variable for staraging edge length
    epsilon = 0.5                                                   
    # Error
    all_edges = []                                                  # Store the length of different edges

    ob = bpy.context.edit_object
    me = ob.data                                                    # Read obj
    bm = bmesh.from_edit_mesh(me)                                   # Read meshs
    bpy.ops.mesh.remove_doubles()                                   # Remove double
    
    # Merger nearby vertices
    for vert1 in bm.verts:
        for vert2 in bm.verts:
            bpy.ops.mesh.select_all(action='DESELECT') 
            if vert1.index != vert2.index and mesh.mesh.cal_length(vert1.co,vert2.co) <= epsilon:
                vert1.select = True
                vert2.select = True 
                bpy.ops.mesh.merge(type='CENTER')
    
    # Record vertices co for deviding groups
    gb_group_co = []
    for face in bm.faces:
        gb_group_co.append([])
        for vert in face.verts:
            gb_group_co[-1].append(vert.co)
    
    # Subdide GBs 
    for edge in bm.edges:
        ini_all_edges.append(round(edge.calc_length(),4))
    ini_all_edges.sort()
    ini_edges = ini_all_edges[1::2]
    edge_counts = Counter(ini_edges)
    all_edges.append(edge_counts.most_common()[0][0])
    vertical_edges = [i for i in ini_edges if i != all_edges[0]]

    if len(vertical_edges) != 0:
        for i in range(len(vertical_edges)):
                all_edges.append(vertical_edges[i])

    for select_edge in all_edges:
        mesh.mesh.cut_edge(me, bm, select_edge, mesh_width, epsilon)
        
    # Devide GB mesh into several groups   
    vert_group_key = np.full(len(bm.verts), False)
    for gb_index in range(len(gb_group_co)):
        bpy.ops.mesh.select_all(action='DESELECT') 
        gb_group_verts = gb_group_co[gb_index]                # four points in a GB plane
        gb_hull = mesh.mesh.cal_hull(gb_group_verts,epsilon)
        group = bpy.context.object.vertex_groups.new()        # creat a new group
        if gb_index < 9:                                      # name the group
            vertex_g_name = ("0" + str(gb_index+1))
        else:
            vertex_g_name = str(gb_index+1)
        group.name = vertex_g_name
        for vert in bm.verts:                                 # check whether the new verts of GB messh belongs to the group
            if mesh.mesh.in_hull(vert.co,gb_hull) and vert_group_key[vert.index] == False:
                vert_group_key[vert.index] = True
                vert.select = True
        bmesh.update_edit_mesh(me, True) 
        bpy.ops.object.vertex_group_assign()
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', ngon_method='BEAUTY')
    bpy.ops.mesh.vertices_smooth()
    return {'FINISHED'}

def tip_hull(self, context):
    print("")
    print("Mesh in tip execution")
                
    # read the mesh file
    ob = bpy.context.edit_object
    me = ob.data                                                        # Read obj
    bm = bmesh.from_edit_mesh(me)

    bpy.ops.mesh.select_all(action='DESELECT')

    # create a vertex list from the objects vertices
    idx = object_idx(bpy.context.scene.gb_mesh_panel_props.target_object)            # Rhett changed 20210522
    target_object_verts = bpy.context.scene.objects[idx].data.vertices
    poslist = []
    targetWorldMatrix = bpy.context.scene.objects[idx].matrix_world
    for vert in target_object_verts:
        poslist.append(targetWorldMatrix @ vert.co)
    
    # generate the hull from pos file
    hull_pts = []
    hull = ConvexHull(poslist)
    hull_indices = hull.vertices
    surf_vert_num = len(hull_indices)
    surf_points_gap = floor(surf_vert_num/100)
    for i in hull_indices:
        if i % surf_points_gap == 0:
            hull_pts.append(poslist[i])
    #print(hull_pts)
    tip_hull = Delaunay(hull_pts)
    
    for vert in bm.verts:
        vert.select = not mesh.mesh.in_hull(vert.co,tip_hull)
                
    bpy.ops.mesh.delete(type='VERT')
    bpy.ops.mesh.select_all(action='SELECT')
        
    return {'FINISHED'}

def Dcom(self, context):
    print("")
    print("DCOM execution")
    wasEdit = False
    if bpy.context.mode == 'EDIT_MESH':
        bpy.ops.object.editmode_toggle()
        wasEdit = True
    dcom_radius = bpy.context.scene.gb_mesh_panel_props.dcom_radius

    idx = object_idx(bpy.context.scene.gb_mesh_panel_props.target_object)            # Rhett changed 20210522
    target_object_verts = bpy.context.scene.objects[idx].data.vertices
    
    vlist = []
    targetWorldMatrix = bpy.context.scene.objects[idx].matrix_world
    
    for vert in target_object_verts:
        vlist.append(targetWorldMatrix @ vert.co)

    #ACTUAL DCOM STEPS 
    worldMatrix = bpy.context.object.matrix_world
    # relaxation along vertex normals
    for vert in bpy.context.object.data.vertices:
        # moving selected vertices along their normal
        if vert.select:
            nor = vert.normal
            ver = worldMatrix @ vert.co
            dcom = Vector((0,0,0))
            count = 0
            for tVert in vlist:
                dVector = tVert - vert.co
            
                if dVector.length <= dcom_radius:
                    dcom = dcom + dVector
                    count += 1
            if count>0:
                dcom = dcom * (1/count)      
                increment = dcom @ nor
                #print(increment)
                vert.co.x = vert.co.x + nor.x*increment
                vert.co.y = vert.co.y + nor.y*increment
                vert.co.z = vert.co.z + nor.z*increment 

    bpy.ops.object.editmode_toggle()         
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.remove_doubles()
    bpy.ops.mesh.beautify_fill()
    bpy.ops.mesh.vertices_smooth()
    bpy.ops.mesh.delete_loose()
    bpy.ops.object.editmode_toggle()         

    return {'FINISHED'}

def object_idx(object_name):
    for index, object in enumerate(bpy.context.scene.objects): #iterate over all objects
        if object.name == object_name:
            return index

def ply_color(self, context):
    obj = context.object
    objmat = bpy.ops.material.new()
    blend.material.set(obj, objmat)
    return {}