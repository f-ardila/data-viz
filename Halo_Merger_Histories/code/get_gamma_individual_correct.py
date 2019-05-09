import time
start_time = time.time()

# Working with merger trees: https://ytree.readthedocs.io/en/latest/Arbor.html

def get_id_given_a(my_tree, scale_factor): 
    
    a_idx = np.where(my_tree["tree", "scale_factor"] == scale_factor)[0]
    
    if len(a_idx)==0:
        return -1
    
    id_at_sfnow = my_tree["tree", "uid"][a_idx]

    return np.ndarray.flatten(id_at_sfnow)[0]

def get_treeidx_given_a(my_tree, scale_factor): 
    
    a_idx = np.where(my_tree["tree", "scale_factor"] == scale_factor)[0]
    
    if len(a_idx)==0:
        return (-1, -1)
    
    id_at_sfnow = my_tree["tree", "uid"][a_idx]
    
    if len(id_at_sfnow) > 0 and halo in id_at_sfnow:
        # if yes then find the tree node now, and then the tree node then 
        treeidx = np.where(halo==id_at_sfnow)[0]
        treenow = (my_tree["tree"][a_idx])[treeidx]
        return (treeidx, treenow[0])
    else:
        return (-1, -1 )
    
def get_tree_given_a(my_tree, scale_factor): 
    a_idx = np.where(my_tree["prog", "scale_factor"] == scale_factor)
    
    if len(a_idx)==0:
        return -1 
    
    if len(my_tree["prog"][a_idx]) > 0:
        return my_tree["prog"][a_idx][0]
    else:
        return -1 
    
def get_alltrees_given_a(my_tree, scale_factor): 
    a_idx = np.where(my_tree["tree", "scale_factor"] == scale_factor)
    
    if len(a_idx)==0:
        return -1 
    
    if len(my_tree["tree"][a_idx]) > 0:
        return my_tree["tree"][a_idx]
    else:
        return -1 
    
import ytree 
import numpy as np
import h5py
from colossus.halo.mass_so import dynamicalTime
from colossus.cosmology import cosmology
import matplotlib.pyplot as plt 
import glob
import pandas as pd
import os 
import sys 

treedir = "/zang/pbehrooz/MDPL2/trees/"

treefile = sys.argv[1]
rank = sys.argv[2]

# if os.path.isfile("/pfs/home/exhakaj/Output/%s.csv" % treefile): continue 

print ("-------- WORKING WITH TREE FILE %s, rank %s ----------- " % (treefile, rank))

a = ytree.load(treedir+treefile)
a.set_selector("max_field_value", "mvir") # the first progenitor descendant will be the most massive one


# ----- OPEN THE HALO CATOALOG YOU WANT TO WORK WITH ------- 
cosmo = cosmology.setCosmology('planck13') # CHANGE THIS LATER 
hlist_dir = "/pfs/home/exhakaj/Hlist_dir/"
sf_now = 0.73330

# ----- GET THE HALO ID'S AND THEIR M200b --------
# halo_now = h5py.File(hlist_dir + "hlist_%.5f.hdf5" % sf_now, "r")["data"]
halo_now = h5py.File(hlist_dir + "hlist_%.5f.hdf5" % sf_now, "r")["data"]


## get dynamical time 
z_now = (1.-sf_now)/sf_now
dyntime = dynamicalTime(z_now, "200m", definition='crossing')*1e9 # in years

# ----- FIND SF 1 TDYN AGO USING SPARTA'S METHOD -------

# ## time now 
time_snapshot = cosmo.age(z_now, inverse = False) *1e9

# ## time back 1tdyn
time_back1dyn = time_snapshot - dyntime

## scale factor itdyn ago
z_back1dyn = cosmo.age(time_back1dyn/1e9, inverse = True) # with inverse = True we compute z(t)
sf_back1dyn = 1/(1+z_back1dyn)


## find the closest snapshot 
catalogs = glob.glob("/zang/pbehrooz/MDPL2/hlists/hlist*")

sf_all = []
for cat in catalogs: 
    sf = cat.split("/")[-1]
    sf = sf.replace("hlist_", "")
    sf = sf.replace(".list", "")
    sf_all.append(float(sf))
    
sf_all = np.array(sf_all)

sf_then = sf_all[np.argmin(abs(sf_all - sf_back1dyn))]
# sf_then = 0.79927

# ----- GET THE ID OF THE HALO 1TDYN AGO --------

ID_now_cat = halo_now["halo_id"]
    

ID_now_tree = []
Idx_now_tree = []

# collect all the halost for the given scale factor in the tree 
for itree, tree in enumerate(a): # a == list of trees for the for loop 
    id_now = get_id_given_a(tree, sf_now) # given_a means given sf 
    
    if id_now == -1: continue 

    ID_now_tree.append(id_now)
    Idx_now_tree.append(itree)

        
        
print ("Finished collecting all of my halos")

Idx_now_tree = np.array(Idx_now_tree)
ID_now_tree = np.array(ID_now_tree)

ID_then = []
ID_now = []
Gamma = []
M200b_now = []
M200b_then = []

for ihalo, halo in enumerate(ID_now_cat): 
    
    # find in which tree this halo lies
    idx_all_tree = np.where(ID_now_tree == halo)[0]
    
#     print (idx_all_tree)
    if not np.isscalar(idx_all_tree):
        if len(idx_all_tree) == 0: 
            continue 
        else: 
            idx_all_tree = idx_all_tree[0] # get the most massive descendant/progenitor to be the main one 
    
    idx_tree = Idx_now_tree[idx_all_tree]
    
#     print (Idx_now_tree[idx_all_tree], ID_now_tree[idx_all_tree])
    
    Treenow = get_alltrees_given_a(a[idx_tree], sf_now)
    
    for treenow in Treenow:
        treethen = get_tree_given_a(treenow, sf_then)
    
        if treethen == -1: continue

        # get the id of the tree node then 
        id_then = treethen["uid"]
        id_now = treenow["uid"]
        # get m200b of the tree node then 
        m200b_then = float(treethen["M200b"])

        # get m200b now 
        m200b_now = float(treenow["M200b"])

        # compute gamma 
        gamma = np.log10(m200b_now/m200b_then)/np.log10(sf_now/sf_then)

        Gamma.append(gamma)
        ID_then.append(id_then)
        ID_now.append(id_now)
        M200b_now.append(m200b_now)
        M200b_then.append(m200b_then)
    
data = {'id_now':ID_now, 'id_then':ID_then, 'gamma_sp': Gamma, "m200b_then":M200b_then, "m200b_now":M200b_now}

df = pd.DataFrame(data=data)

df.to_csv("/pfs/home/exhakaj/Output_m200b/%s.csv" % treefile, sep=',', index = None)

print ("Done with this part! %s" % treefile)        

print ("My program took %s to run" % (time.time() - start_time))
