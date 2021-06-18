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
import os

# bpy types used
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, FloatProperty, EnumProperty
from bpy.types import Operator

# Own pkgs
from . import GBMesh_opexe as opexec

HALO_IMG_PATH = os.path.dirname(__file__)+"/atomtex.png"

# === Operator classes ===
class VIEW3D_OT_pospath_button(Operator, ImportHelper):
    """Select POS file from dialogue"""
    bl_idname = "gbmesh.import_pospath"
    bl_label = "Select .pos file"

    # ImportHelper mixin class uses this
    filename_ext = ".pos"

    filter_glob = StringProperty(
            default="*.pos",
            options={'HIDDEN'},
            )

    def execute(self, context):
        props = context.scene.pos_panel_props
        # set pos filename
        props.pos_filename = self.filepath
        return {'FINISHED'}

class VIEW3D_OT_rngpath_button(Operator, ImportHelper):
    """Select RNG file from dialogue"""
    bl_idname = "gbmesh.import_rngpath"
    bl_label = "Select .rng file"

    # ImportHelper mixin class uses this
    filename_ext = ".rng"

    filter_glob = StringProperty(
            default="*.rng",
            options={'HIDDEN'},
            )

    def execute(self, context):
        props = context.scene.pos_panel_props
        # set rng filename
        props.rng_filename = self.filepath
        return {'FINISHED'}

# class VIEW3D_OT_load_posrng_button(Operator):
#     """Read and load POS/RNG files into memory"""
#     bl_idname = "gbmesh.load_posrng"
#     bl_label = "Load POS/RNG files"

#     def execute(self, context):
#         return opexec.load_posrng(self, context)

class VIEW3D_OT_bake_button(Operator):
    """Bake POS data to object"""
    bl_idname = "gbmesh.bake_button"
    bl_label = "Bake to object"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        return (area == 'VIEW_3D') and (mode == 'OBJECT')

    def execute(self, context):
        return opexec.bake(self, context)

class VIEW3D_OT_clear_button(Operator):
    """Clears all meshes in the scene"""
    bl_idname = "gbmesh.clear_button"
    bl_label = "Clear all objects"

    def execute(self, context):
        return opexec.clear(self, context)

class VIEW3D_OT_import_obj_button(Operator, ImportHelper):
    """Select OBJ file from dialogue"""
    bl_idname = "gbmesh.import_obj_button"
    bl_label = "Select .obj file"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        return (area == 'VIEW_3D') and (mode == 'OBJECT')

    # ImportHelper mixin class uses this
    filename_ext = ".obj"

    filter_glob = StringProperty(
            default="*.obj",
            options={'HIDDEN'},
            )

    def execute(self, context):
        props = context.scene.obj_panel_props
        # set import obj filename
        props.import_obj_filename = self.filepath
        return opexec.load_obj(self, context)

class VIEW3D_OT_export_obj_button(Operator):
    """Select OBJ file from dialogue"""
    bl_idname = "gbmesh.export_obj_button"
    bl_label = "Save .obj file"

    def execute(self, context):
        return opexec.save_obj(self, context)

class VIEW3D_OT_import_ply_button(Operator, ImportHelper):
    """Select PLY file from dialogue"""
    bl_idname = "gbmesh.import_ply_button"
    bl_label = "Select .ply file"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        return (area == 'VIEW_3D') and (mode == 'OBJECT')

    # ImportHelper mixin class uses this
    filename_ext = ".ply"

    filter_glob = StringProperty(
            default="*.ply",
            options={'HIDDEN'},
            )

    def execute(self, context):
        props = context.scene.obj_panel_props
        # set import ply filename
        props.import_ply_filename = self.filepath
        return opexec.load_ply(self, context)

class VIEW3D_OT_remove_duplivert(Operator):
    """Remove vertex duplication object on active object"""
    bl_idname = "gbmesh.remove_duplivert"
    bl_label = "Remove vertex object"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        obj = context.object
        return (area == 'VIEW_3D') and (mode == 'OBJECT') and (obj is not None) and (obj.vistype == 'DUPLI')

    def execute(self, context):
        return opexec.remove_duplivert(self, context)

class VIEW3D_OT_dupli_vert(Operator):
    """Duplicate icosphere on vertices of currently active object"""
    bl_idname = "gbmesh.add_duplivert"
    bl_label = "Add vertex object"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        obj = context.object
        return (area == 'VIEW_3D') and (mode == 'OBJECT') and (obj is not None)

    def execute(self, context):
        return opexec.dupli_vert(self, context)

class VIEW3D_OT_add_halo_material(Operator):
    """Apply halo material to currently selected object"""
    bl_idname = "gbmesh.add_halomat"
    bl_label = "Add halo material"

    # path to billboard texture
    halo_img_path = StringProperty(
            description = "Image to use for halo texture",
            default = HALO_IMG_PATH
            )

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        obj = context.object
        return (area == 'VIEW_3D') and (mode == 'OBJECT') and (obj is not None)

    def execute(self, context):
        return opexec.add_halo_material(self, context)

class VIEW3D_OT_remove_halo_material(Operator):
    """Remove halo material on currently selected object"""
    bl_idname = "gbmesh.remove_halomat"
    bl_label = "Remove halo material"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        obj = context.object
        return (area == 'VIEW_3D') and (mode == 'OBJECT') and (obj is not None) and (obj.vistype == 'HALO')

    def execute(self, context):
        return opexec.remove_halo_material(self, context)

class VIEW3D_OT_add_bounding_box(Operator):
    """Add a bounding box to current data"""
    bl_idname = "gbmesh.add_bound_box"
    bl_label = "Add bounding box"

    # Boundbox padding
    padding = FloatProperty(
            name="Padding",
            description="Bounding box padding",
            default=0.5,
            min=0.0, max=100.0,
            )

    def execute(self, context):
        return opexec.add_bounding_box(self, context)

class VIEW3D_OT_pointcloud_add(Operator):
    """Create pointcloud visualisation for the active object"""
    bl_idname = "gbmesh.pointcloud_add"
    bl_label = "Create pointcloud"

    @classmethod
    def poll(cls, context):
        area = context.area.type
        mode = context.mode
        objtype = None
        if context.object is not None:
            objtype = context.object.type
        return (area == 'VIEW_3D') and (mode == 'OBJECT') and (objtype == 'MESH')

    def execute(self, context):
        return opexec.pointcloud_add(self, context)

class VIEW3D_OT_scale_child(Operator):
    """Scale child of currently selected object"""
    bl_idname = "gbmesh.scale_child"
    bl_label = "Scale vertex object"

    @classmethod
    def poll(cls, context):
        obj = context.object
        return (obj is not None) and obj.children

    def execute(self, context):
        return opexec.scale_child(self, context)


## ++++ The Following is for GBMesh ++++
class OBJECT_OT_quad_mesh(Operator):
    bl_idname = "gbmesh.quad_mesh_button"
    bl_label = "Quad mesh"
    
    def execute(self,context):
        return opexec.quad_mesh(self, context)
        
### actual GB in tip full operation
class OBJECT_OT_tip_hull(Operator):
    bl_idname = "gbmesh.tip_hull_button"
    bl_label = "Mesh in tip"
    
    def execute(self,context):
        return opexec.tip_hull(self, context)

### actual DCOM operation                                                     Copy from Peter Johann Felfer
class OBJECT_OT_dcom(Operator):
    bl_idname = "gbmesh.dcom_button"
    bl_label = "DCOM"
    
    def execute(self,context):
        return opexec.Dcom(self, context)
