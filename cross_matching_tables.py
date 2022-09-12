#!/usr/bin/env python
# coding: utf-8

# # Matching astrometrical catalogs

# In[36]:


"""import libraries"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


import glob


# In[37]:


"""define tables path and names"""

dirpath = './files/'

files = glob.glob(dirpath + '*.csv')


"""columns use for matching"""

c2m = ['RAJ2000','DEJ2000']


"""Threshold"""

th = 10./3600.


# In[38]:


"""Reading and sorting tables"""

df = [pd.read_csv(f) for f in glob.glob(dirpath + '*.csv')]

"""Sorting"""

df.sort(key=len, reverse=False)

coords1, coords2 = df[0][c2m], df[1][c2m]
coords1.columns = coords2.columns = ['x', 'y']

coords1 = coords1.sort_values(by=['x'])
coords2 = coords2.sort_values(by=['x'])


# $$d = \sqrt{\left[ (\alpha_{A} - \alpha_{B})\cdot\cos(\delta_{B}) \right] ^2 + (\delta_{A} - \delta_{B})^2}$$

# In[39]:


def dst(t1,t2):
    return np.sqrt(((t1.x - t2.x)*np.cos(np.pi*t2.y/180.0))**2 + (t1.y-t2.y)**2)


# In[40]:


get_ipython().run_cell_magic('time', '', '\n"""construct the distance matrix"""\n\ndst_m = coords1.apply(dst, t2=coords2, axis=1)')


# In[41]:

"""matrix plot"""

fig = plt.figure(figsize=(10,10))
plt.imshow(dst_m, cmap='hot', interpolation='None')
fig.savefig("dtsm-imshowAB.png", dpi=100)


# In[42]:


get_ipython().run_cell_magic('time', '', '"""threshold cut"""\n\ndst_nm = dst_m.values\n\nmd_m = dst_m[(dst_m < th)]')


# In[43]:


fig = plt.figure(figsize=(10,10))
plt.imshow(md_m, cmap='hot', aspect='equal', interpolation='None')
fig.savefig("dtsm-imshowABthr.png", dpi=100)


# In[44]:


ra1, de1 = df[0][c2m[0]],df[0][c2m[1]]
ra2, de2 = df[1][c2m[0]],df[1][c2m[1]]


# In[49]:


get_ipython().run_cell_magic('time', '', "\ncl1 = ['nmA']*len(ra1)\ncl2 = ['nmB']*len(ra2)\n\n\nfor i, ra, de in zip(range(len(ra1)),ra1,de1):\n    d = np.sqrt(((ra-ra2)*np.cos(de*np.pi/180.))**2.+(de-de2)**2.)\n\n    if(len(np.where(d < th)[0]) > 0):\n        k  = np.argwhere((d < th) & (d == d.min()))[0][0]\n        dk = np.sqrt(((ra1-ra2[k])*np.cos(de2[k]*np.pi/180.))**2.+(de1-de2[k])**2.)\n        el = np.argwhere((dk<th) & (dk == dk.min()))[0][0]\n        cl1[el] = i\n        cl2[k]  = i\n        \ndf[0]['idMA'] = cl1\ndf[1]['idMB'] = cl2")


# In[46]:


pdm = pd.merge(df[0],df[1], how = 'outer', left_on=['idMA'], right_on=['idMB'])

pdm[c2m[0]] = pdm[[c2m[0]+'_x',c2m[0]+'_y']].mean(axis=1)
pdm[c2m[1]] = pdm[[c2m[1]+'_x',c2m[1]+'_y']].mean(axis=1)


# In[47]:


pdm.to_csv('match-out-catalog.csv', index = True, index_label = 'id')


# In[48]:


corr_ra=pdm[c2m[0]+'_x'].corr(pdm[c2m[0]+'_y'])
corr_de=pdm[c2m[1]+'_x'].corr(pdm[c2m[1]+'_y'])
corr_rajy=pdm[c2m[0]].corr(pdm[c2m[0]+'_y'])
corr_dejy=pdm[c2m[1]].corr(pdm[c2m[1]+'_y'])
corr_rajx=pdm[c2m[0]].corr(pdm[c2m[0]+'_x'])
corr_dejx=pdm[c2m[1]].corr(pdm[c2m[1]+'_x'])

corr_ra, corr_de, corr_rajy, corr_dejy, corr_rajx, corr_dejx


# In[ ]:




