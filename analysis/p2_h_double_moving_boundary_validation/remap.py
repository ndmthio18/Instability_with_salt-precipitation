"""Conservative re-representation; never converts positive crystal by a cutoff."""
import numpy as np

def overlap_remap(old_faces, inventory, new_faces):
    old_faces=np.asarray(old_faces,float);new_faces=np.asarray(new_faces,float)
    inventory=np.asarray(inventory,float)
    if not (np.all(np.diff(old_faces)>0) and np.all(np.diff(new_faces)>0)):
        raise ValueError('nonpositive cell volume')
    if abs(old_faces[0]-new_faces[0])+abs(old_faces[-1]-new_faces[-1])>1e-14:
        raise ValueError('remap domains differ')
    overlap=np.maximum(0,np.minimum(new_faces[1:,None],old_faces[None,1:])
                         -np.maximum(new_faces[:-1,None],old_faces[None,:-1]))
    return overlap @ (inventory/np.diff(old_faces))

def convert_connected_zero_prefix(liquid_volume,liquid_salt,crust_volume,solid_volume,crust_salt):
    v=np.asarray(crust_volume,float);s=np.asarray(solid_volume,float);q=np.asarray(crust_salt,float)
    if np.any(s<0) or np.any(s>v):raise ValueError('invalid crystal inventory')
    count=0
    while count<len(s) and s[count]==0:count+=1
    # Exact conservative transfer of the liquid already present, not a reaction.
    return (np.r_[liquid_volume,v[:count]],np.r_[liquid_salt,q[:count]],
            v[count:].copy(),s[count:].copy(),q[count:].copy(),count)
