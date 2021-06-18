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

from bpy.types import PropertyGroup, Panel
from bpy.props import StringProperty, EnumProperty, FloatProperty

from . import GBMesh_opexe as opexec

DEFAULT_COLOR = (0, 0.144, 0.554)

# === GBMesh-specific object RNA properties ===
# Define GBMesh-specific RNA props for every object
# Set to True for top level objects (eg elements)

# Defines object type in GBMesh framework
# Default: BLENDER for objects independent of GBMesh
dtypes = [('BLENDER', "Blender",  "Blender"),
          ('DATA',    "Dataset",  "Dataset"),
          ('BOUND',   "Boundbox", "Boundbox")]
bpy.types.Object.datatype = EnumProperty(
        name = "Type of object (GBMesh)",
        items = dtypes,
        default = 'BLENDER'
        )

# Type of visualisation applied
vtypes = [('NONE',  "None",  "None"),
          ('HALO',  "Halo",  "Halo"),
          ('DUPLI', "Dupli", "Dupli")]
bpy.types.Object.vistype = EnumProperty(
        name = "Type of visualisation (GBMesh)",
        items = vtypes,
        default = 'NONE'
        )

# === General panel properties ===
class VIEW3D_PT_pos_panel_props(PropertyGroup):
    """ POS reader panel property group

    Properties:
    pos_filename -- POS file path
    rng_filename -- RNG file path
    plot_type -- Enumerator in ['ISO', 'EA', 'ION']
                 Plot by isotope, atom, or ion
    atoms, rngs, ions -- Enumerators for atoms/rngs/ions loaded from files
    """

    pos_filename = StringProperty(
            name = "",
            description = "Input .pos file",
            default = ""
        )

    rng_filename = StringProperty(
            name = "",
            description = "Input .rng file",
            default = ""
        )

    plot_options = [('EA', "Atomic", "Atomic"), ('ION', "Ionic", "Ionic"), ('ISO', "Isotopic", "Isotopic")]
    plot_type = EnumProperty(name="Bake options", items=plot_options)

class VIEW3D_PT_obj_panel_props(PropertyGroup):
    """ OBJ reader panel property group

    Properties:
    import_obj_filename -- OBJ file path for importing
    export_obj_filename -- OBJ file path for exporting
    import_ply_filename -- ply file path for importing
    """

    import_obj_filename = StringProperty(
            name = "",
            description = "Import .obj file",
            default = ""
        )
    
    # import_obj_filepath = StringProperty(
    #         name = "",
    #         description = "Folder for imported .obj file",
    #         default = ""
    #     )

    export_obj_filename = StringProperty(
            name = "",
            description = "Export .obj file",
            default = ""
        )

    import_ply_filename = StringProperty(
            name = "",
            description = "Import .ply file",
            default = ""
        )


class VIEW3D_PT_gb_mesh_panel_props(PropertyGroup):
    """ GB Mesh panel property group

    Properties:
    quad_mesh_edge_length -- Minimum edge length of a quad mesh
    target_object -- The target object that is used for the grain boundary meshing
    dcom_radius -- Radius around a vertex within which the center of mass is calculated         Copy from Peter Johann Felfer
    """

    quad_mesh_edge_length = FloatProperty(
            name = "Mesh length", default = 3, 
            min = 1, max = 20)
 
    target_object = StringProperty( 
                name="gbmesh target", 
                description="GB mesh target object", 
                default= ""
                )

    dcom_radius = FloatProperty(
            name = "DCOM radius", default = 5, 
            min = 0, max = 100)

class VIEW3D_PT_import_panel(Panel):
# ++++ Creates a Panel in the scene context of the properties editor ++++
    bl_idname = "Import_PT_Panel"
    bl_label = "Import APT data"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "GBMesh"
    
    def draw(self, context):
        layout = self.layout
        props = context.scene.pos_panel_props

        col = layout.column(align=True)
        col.operator("GBMesh.import_pospath")
        col.operator("GBMesh.import_rngpath")
        # col.operator("GBMesh.load_posrng")

        col = layout.column(align=True)
        col.prop(props, "plot_type", text="")
        col.operator("GBMesh.bake_button")

        row = layout.row()
        row.operator("GBMesh.clear_button")
        
    # Only display in object mode
    @classmethod
    def poll(cls, context):
        mode = context.mode
        return mode == 'OBJECT'

class VIEW3D_PT_obj_panel(Panel):
# ++++ Creates a Panel in the scene context of the properties editor ++++
    bl_idname = "Object_PT_Panel"
    bl_label = "Import and export mesh data"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "GBMesh"
    
    def draw(self, context):
        layout = self.layout
        props = context.scene.obj_panel_props

        col = layout.column(align=True)
        col.operator("GBMesh.import_obj_button")
        col.operator("GBMesh.export_obj_button")
        #col.operator("GBMesh.import_ply_button")
                
    # Only display in object mode
    @classmethod
    def poll(cls, context):
        mode = context.mode
        return mode == 'OBJECT'


class VIEW3D_PT_data_visualisation(Panel):
    """Visualisation panel viewed in object mode with an AtomBlend-generated dataset selected"""
    bl_idname = "Visualisation_PT_Panel"
    bl_label = "Visualisation"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "GBMesh"

    def draw(self, context):
        layout = self.layout
        props = context.scene.obj_panel_props

        col = layout.column(align=True)
        col.operator("GBMesh.import_ply_button")

        obj = context.object

        if obj is None or obj.datatype != 'DATA':
            col.label(text="[no dataset selected]")
        elif obj.vistype == 'HALO' and has_halo(obj):
            #--- halo material ---
            # active object already has applied halo material, show edit props
            col.label(text="Edit halo material:")
            mat = context.object.active_material
            halo = mat.halo

            subrow = col.row(align=True)
            subrow.prop(mat, "diffuse_color", text="")
            subrow.prop(halo, "size")
            col.operator("GBMesh.remove_halomat")
        elif obj.vistype == 'DUPLI' and has_duplivert(obj):
            #--- duplivert ---
            # active object has duplivert applied, show edit props
            dupli = obj.children[0]
            mat = dupli.active_material

            col.label(text="Edit vertex object:")
            subrow = col.row(align=True)
            subrow.prop(mat, "diffuse_color", text="")
            subrow.operator("GBMesh.scale_child")
            col.operator("GBMesh.remove_duplivert")
        elif obj.datatype == 'DATA':
            # no visualisation applied yet
            #col.operator("GBMesh.add_halomat")
            col.operator("GBMesh.add_duplivert")

        # === Boundbox ===
        row = layout.row()
        col = layout.column(align=True)

        if (obj is not None) and obj.datatype == 'BOUND' and is_bound(obj):
            # boundbox selected
            mod = obj.modifiers[0]
            mat = obj.active_material
            # TODO add remove boundbox button
            col.label(text="Edit bounding box:")
            col.prop(mod, "thickness", text="Thickness")
            col.prop(mat, "diffuse_color", text="")
        else:
            boundbox_props = col.operator("GBMesh.add_bound_box")
            col.prop(boundbox_props, "padding")
            #bpy.ops.atomblend.add_bound_box.padding

# === Helper functions ===
def has_halo(obj):
    mat = obj.active_material
    return (mat is not None) and (mat.type == 'HALO')

def has_duplivert(obj):
    if obj is None:
        return False
    return obj.children and (obj.instance_type == 'VERTS')

def is_bound(obj):
    if obj is None:
        return False
    mat = obj.active_material
    mods = obj.modifiers
    return mods and (mat is not None)

# ++++ The following is for GBMesh ++++
### creating the UI panel
class VIEW3D_PT_gb_mesh_panel(Panel):
    bl_idname = "GB_mesh_PT_Panel"
    bl_label = "Process GB"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "GBMesh"
    bl_context = "mesh_edit"
    
    @classmethod
    def poll(self,context):
        if context.object and context.object.type == "MESH":
            return True
    
    def draw(self, context):
        layout = self.layout
        props = context.scene.gb_mesh_panel_props

        col = layout.column(align=True)
        col.prop(props,"quad_mesh_edge_length", text = "Mesh length (nm)") #minimum edge length of a quad mesh
        col.operator("GBMesh.quad_mesh_button", text = "Quad mesh") #execute quad mesh

        col = layout.column(align=True)
        col.prop(props,"target_object", text="Target") #draw dropdown box on panel
        col.operator("GBMesh.tip_hull_button", text = "Mesh in tip") #execute filter mesh by tip hull

        col = layout.column(align=True)
        col.prop(props,"dcom_radius", text = "DCOM radius (nm)")
        col.operator("GBMesh.dcom_button", text = "DCOM")   #execute dcom to the mesh 