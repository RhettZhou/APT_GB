# =============================================================================
# (C) Copyright 2021
# Max-Planck-Institut für Eisenforschung GmbH
# The University of Alabama
# The University of Sydney
# =============================================================================
# File:   __init__.py
# Date:   2021-05-20
# Author: Xuyang (Rhett) Zhou
#         Adapted from AtomBlend by Varvara Efremova (The University of Sydney)
#
# Description:
# GBMesh addon initialisation and UI definition

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions are met:
    
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# =============================================================================

bl_info = {
    "name": "GBMesh",
    "author": "Rhett Xuyang Zhou", "Varvara Efremova"
    "version": (0, 1),
    "blender": (2, 80, 0),
    "location": "View3D > UI",
    "description": "Atom Probe data GB mesh plugin",
    "warning": "",
    "doc_url": "",
    "category": "3D View",
}

# ++++ Addon reload support ++++
import bpy

from . GBMesh_op    import VIEW3D_OT_pospath_button, VIEW3D_OT_rngpath_button, VIEW3D_OT_bake_button,\
                            VIEW3D_OT_clear_button, VIEW3D_OT_import_obj_button, VIEW3D_OT_export_obj_button, VIEW3D_OT_import_ply_button,\
                            VIEW3D_OT_remove_duplivert,  VIEW3D_OT_dupli_vert, VIEW3D_OT_add_halo_material,\
                            VIEW3D_OT_remove_halo_material, VIEW3D_OT_add_bounding_box, VIEW3D_OT_pointcloud_add, VIEW3D_OT_scale_child,\
                            OBJECT_OT_quad_mesh, OBJECT_OT_tip_hull, OBJECT_OT_dcom

from . GBMesh_panel import VIEW3D_PT_pos_panel_props, VIEW3D_PT_obj_panel_props, VIEW3D_PT_gb_mesh_panel_props,\
                            VIEW3D_PT_import_panel, VIEW3D_PT_obj_panel, VIEW3D_PT_gb_mesh_panel,VIEW3D_PT_data_visualisation

classes = (
    VIEW3D_OT_import_obj_button,
    VIEW3D_OT_export_obj_button,
    VIEW3D_OT_import_ply_button,
    VIEW3D_OT_scale_child,
    VIEW3D_OT_pointcloud_add,
    VIEW3D_OT_add_bounding_box,
    VIEW3D_OT_add_halo_material,
    VIEW3D_OT_remove_halo_material,
    VIEW3D_OT_dupli_vert,
    VIEW3D_OT_remove_duplivert,
    VIEW3D_OT_clear_button,
    VIEW3D_OT_bake_button,
#    VIEW3D_OT_load_posrng_button,
    VIEW3D_OT_rngpath_button,
    VIEW3D_OT_pospath_button,
    VIEW3D_PT_pos_panel_props,
    VIEW3D_PT_obj_panel_props,
    VIEW3D_PT_gb_mesh_panel_props,
    VIEW3D_PT_import_panel,
    VIEW3D_PT_obj_panel,
    VIEW3D_PT_gb_mesh_panel,
    VIEW3D_PT_data_visualisation,
    OBJECT_OT_dcom,
    OBJECT_OT_tip_hull,
    OBJECT_OT_quad_mesh,
)

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    bpy.types.Scene.pos_panel_props = bpy.props.PointerProperty(type=VIEW3D_PT_pos_panel_props)
    bpy.types.Scene.obj_panel_props = bpy.props.PointerProperty(type=VIEW3D_PT_obj_panel_props)
    bpy.types.Scene.gb_mesh_panel_props = bpy.props.PointerProperty(type=VIEW3D_PT_gb_mesh_panel_props)

def unregister():
    from bpy.utils import unregister_class
    del bpy.types.Scene.pos_panel_props
    del bpy.types.Scene.obj_panel_props
    del bpy.types.Scene.gb_mesh_panel_props
    for cls in reversed(classes):
        unregister_class(cls)

if __name__ == "__main__":
    register()