'''
Created on Oct 16, 2016

'''


import os
import time

import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from matplotlib import cm
from sklearn import manifold

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import numpy as np


class data_brower:
    
    def __init__(self, x_coords, y_coords, z_coords, data_name_list):
        self.X=np.array(x_coords)
        self.Y=np.array(y_coords)
        self.Z=np.array(z_coords)
        self.data_names=data_name_list
        self.fig = None
        self.ax = None
        self.be_zoom=False

     
        
    def get_3dcoord(self, x, y):
        s = self.ax.format_coord(x,y)
        print(s)
        x = ''
        y = ''
        z = ''
        for i in range(s.find('x')+2,s.find('y')-2):
            x = x+s[i]
        for i in range(s.find('y')+2,s.find('z')-2):
            y = y+s[i]
        for i in range(s.find('z')+2,len(s)):
            z = z+s[i]
        return float(x),float(y),float(z)
    
    
    def zoom(self,event):
        base_scale = 1.1
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        cur_zlim = self.ax.get_zlim()


        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1 / base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
#             print event.button

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
        new_z = (cur_zlim[1] - cur_zlim[0]) * scale_factor
        x_coord,y_coord,z_coord=self.get_3dcoord(event.xdata,event.ydata)
        
        relx = (cur_xlim[1] - x_coord)/(cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - y_coord)/(cur_ylim[1] - cur_ylim[0])
        relz = (cur_zlim[1] - z_coord)/(cur_zlim[1] - cur_zlim[0])

        self.ax.set_xlim([x_coord - new_width * (1-relx), x_coord + new_width * (relx)])
        self.ax.set_ylim([y_coord - new_height * (1-rely), y_coord + new_height * (rely)])
        self.ax.set_zlim([z_coord - new_z * (1-relz), z_coord + new_z * (relz)])
        self.be_zoom=True
#         self.ax.figure.canvas.draw()
        self.update_position(event)

    def onPress(self,event):
        self.be_zoom=False
        
    def onpick(self, event):
        if self.be_zoom==True:
            return
        print("onpick")
        N = len(event.ind)
        if not N:
            return True
#         print event.mouseevent
#         print event.artist
#         print event.ind
        N = len(event.ind)
        if not N:
            return True
        x_coord,y_coord,z_coord=self.get_3dcoord(event.mouseevent.xdata,event.mouseevent.ydata)
        distances = np.hypot(x_coord - self.X[event.ind], y_coord - self.Y[event.ind],z_coord-self.Z[event.ind])
#         print distances
        
        indmin = distances.argmin()
        dataind = event.ind[indmin]

        self.lastind = dataind
        self.update_point()
        
        
    def display(self,protein_name):
        self.fig = plt.figure('decoys brower')
        self.ax = self.fig.gca(projection='3d')
        
        self.ax.set_title(protein_name+' folding landscape display')
        triang = mtri.Triangulation(self.X, self.Y)
        self.ax.plot_trisurf(triang, self.Z, cmap=cm.jet, linewidth=0.2)
#         self.text = self.ax.text(0.05, 0.95, 0.05,'selected: none', transform=self.ax.transAxes, va='top')
        
        self.ax.scatter(self.X, self.Y, self.Z, s=5, alpha=0.0, visible=True, picker=5)       
        
        self.selected_inds = []
        self.selected_sc=self.ax.scatter([self.X[0]], [self.Y[0]], [self.Z[0]], s=20, alpha=0.4,c='yellow', visible=True)
        self.selected_inds.append(0)
        
        self.selected_labels = []
        x2, y2, _ = proj3d.proj_transform(self.X[0], self.Y[0], self.Z[0], self.ax.get_proj())
        label=plt.annotate(
            self.data_names[0], 
            xy = (x2, y2), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',size='small',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        self.selected_labels.append(label)
        
#         self.fig.canvas.mpl_connect('button_release_event',self.onpick)
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('motion_notify_event', self.update_position)
        self.fig.canvas.mpl_connect('scroll_event', self.zoom)
        self.fig.canvas.mpl_connect('button_press_event',self.onPress)
        plt.show()
        time.sleep(60)
        plt.savefig(protein_name+'.png')
        plt.close()
        
        
    def update_position(self,event):
        
        label_X, label_Y, _ = proj3d.proj_transform(self.X[self.selected_inds], self.Y[self.selected_inds], self.Z[self.selected_inds], self.ax.get_proj())
        for i in range(len(self.selected_inds)):
            label=self.selected_labels[i]
            label.xy = label_X[i],label_Y[i]
            label.update_positions(self.fig.canvas.renderer)
        self.fig.canvas.draw()  
         
         
    def update_point(self):
        if self.lastind is None:
            return
        
        dataind = self.lastind
#         print dataind
        not_find=True
        for i in range(len(self.selected_inds)):
            if self.selected_inds[i]==dataind:
                self.selected_inds.remove(self.selected_inds[i])
                self.selected_labels[i].remove()
                self.selected_labels.remove(self.selected_labels[i])
                not_find=False
                break
        if not_find:
            self.selected_inds.append(dataind)
            
            x2, y2, _ = proj3d.proj_transform(self.X[dataind], self.Y[dataind], self.Z[dataind], self.ax.get_proj())
            label=plt.annotate(
                self.data_names[dataind], 
                xy = (x2, y2), xytext = (-15, 15),
                textcoords = 'offset points', ha = 'right', va = 'bottom', size='small',
                bbox = dict(boxstyle = 'Round,pad=0.5', fc = 'yellow', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            self.selected_labels.append(label)

        paths=np.column_stack((self.X[self.selected_inds], self.Y[self.selected_inds], self.Z[self.selected_inds]))
#         print paths
        self.selected_sc.remove()
        self.selected_sc=self.ax.scatter(self.X[self.selected_inds], self.Y[self.selected_inds], self.Z[self.selected_inds], s=20, alpha=0.4,
                                c='yellow', visible=True)
        self.fig.canvas.draw()
        
        
        
    
def RMSD(model_a,model_b,be_all):
    super_imposer = Bio.PDB.Superimposer()
    ref_atoms = []
    alt_atoms = [] 

    for ref_atom,alt_atom in zip(model_a.get_atoms(),model_b.get_atoms()):
        ref_atoms.append(ref_atom)
        alt_atoms.append(alt_atom)
    if  model_a.id!=model_b.id:
        print(model_a.id,model_b.id)
    assert model_a.id==model_b.id
 
    super_imposer.set_atoms(ref_atoms, alt_atoms)
    super_imposer.apply(alt_atoms)
    return  super_imposer.rms

def loadPdb(native_pdb_name,pdb_path):
    
    p = PDBParser(QUIET = True)
    pdb_path_list=os.listdir(pdb_path)
    pdb_path_list.sort()
    
    pdb_structure_list=[]
    pdb_name_list=[]
    for index in range(len(pdb_path_list)):
        local_pdb_name=pdb_path_list[index]
#         print local_pdb_name.lower()
        if '.pdb' in local_pdb_name.lower():
            
            pstructure=p.get_structure(local_pdb_name, pdb_path+pdb_path_list[index])[0]
            if local_pdb_name == native_pdb_name:
                pdb_name_list.insert(0,local_pdb_name)
                pdb_structure_list.insert(0,pstructure)
            else:
                pdb_name_list.append(local_pdb_name)
                pdb_structure_list.append(pstructure)
    print('load pdb finish.')
    return pdb_structure_list,pdb_name_list  
      
def loadScore(pdb_name,configure_file):
    
    score_dict={}
    print(configure_file)
    score_reader=open(configure_file)
    line=score_reader.readline()
    
    while line:
        line_list=line.split()

        pdb_key=line_list[0].split('/')[-1].split(':')[0]
        score_dict[pdb_key]=float(line_list[1])

        line=score_reader.readline()

    score_reader.close()
    print('load scores finish.')
    return score_dict

def folding_landscape_3d_plot(protein_name,native_pdb_name,pdb_path,configure_file,RMSD_matrix_file_name):
    
    if pdb_path[-1]!='/':
        print('pdb_path')
        pdb_path=pdb_path+'/'
    x_coords = []
    y_coords = []
    z_coords = []
    pdb_label_list=[]
    score_dict=loadScore(protein_name, configure_file)
    [pdb_structure_list,pdb_name_list]=loadPdb(native_pdb_name,pdb_path)
    need_load_flag=True
    if os.path.exists(protein_name+'_coords.cfg'):

        print ('1')
        need_load_flag=False
        pdb_coords_file_reader=open(protein_name+'_coords.cfg')
        line=pdb_coords_file_reader.readline()
        pdb_name_list2=[]
        while line:
            line_list=line.split()
            if len(line_list)>=4:
#                 print('2')
                pdb_name_list2.append(line_list[0])
                x_coords.append(float(line_list[1]))
                y_coords.append(float(line_list[2]))
                z_coords.append(float(line_list[3]))
                line=pdb_coords_file_reader.readline()
                pdb_label_list.append(line_list[0]+' '+line_list[4])
            else:
                need_load_flag=True
                print('3')
                break
        if need_load_flag != True:
            print('4')
            for pdb_name,z_coord in zip(pdb_name_list,z_coords):
                if   pdb_name not in score_dict:
                    print('5',pdb_name)
                    need_load_flag=True
                    break
                elif score_dict[pdb_name]!=z_coord:
                    print ('6',pdb_name,score_dict[pdb_name],z_coord)
                    need_load_flag=True
                    break
    if not need_load_flag:
        print('file have load')
    else:       
        if pdb_structure_list==None and pdb_name_list ==None:
            return
        pdb_label_list=[]
        print('calcaulate  pdb_RMSD start.')
        pdb_RMSD=[] 
        for indexi in range(len(pdb_structure_list)): 
            b=[]
            print(indexi,len(pdb_structure_list))
            for indexj in range(len(pdb_structure_list)): 
                if indexi == indexj:
                    b.append(0.0);
                elif indexi>indexj:
                    b.append(pdb_RMSD[indexj][indexi])
                else:
                    localRMSD=RMSD(pdb_structure_list[indexi],pdb_structure_list[indexj],1)
                    b.append(round(localRMSD,4))
            pdb_RMSD.append(b)
        print('calcaulate  pdb_RMSD finish.')
        
        print('the pdb label combination start.')
        for index in range(len(pdb_RMSD[0])):
            pdb_label_list.append(pdb_name_list[index]+' '+str(pdb_RMSD[0][index]))
        print('the pdb label combination finish.')
        
        print('calcaulate  RMSD_matrix start.')
        RMSD_matrix_writer=open(RMSD_matrix_file_name,'w') 
        RMSD_matrix_writer.write(' '.join(pdb_name_list)+'\n')
        for i in range(len(pdb_RMSD)):
            print(i,len(pdb_RMSD))
            RMSD_matrix_writer.write(pdb_name_list[i])
            for item in pdb_RMSD[i]:
                RMSD_matrix_writer.write(' '+str(item))
            RMSD_matrix_writer.write('\n')
        RMSD_matrix_writer.close()
        print('calcaulate  RMSD_matrix finish.')
        
        # obtain relative coordinate based on protein distance
        print('calcaulate  relative coordinate start.')
        adist = np.array(pdb_RMSD)
        amax = np.amax(adist)
        adist /= amax
        mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        results = mds.fit(adist)
        coords = results.embedding_
        print ('calcaulate  relative coordinate finish.')
        

        x_coords = coords[:,0]
        y_coords = coords[:,1]
        z_coords = []
        

        for index in range(len(pdb_name_list)):
            z_coords.append(score_dict[pdb_name_list[index]])
        print(x_coords)
        print(z_coords)
        print(pdb_name_list)
        pdb_coords_file_writer=open(protein_name+'_coords.cfg','w')
        for pdb_name,x_coord,y_coord,z_coord,rmsd in zip(pdb_name_list,x_coords,y_coords,z_coords,pdb_RMSD[0]):
            pdb_coords_file_writer.write(pdb_name+'\t'+str(x_coord)+'\t'+str(y_coord)+'\t'+str(z_coord)+'\t'+str(rmsd)+'\n')
            
        pdb_coords_file_writer.close()    
    
    brower=data_brower(x_coords,y_coords,z_coords,pdb_label_list)
    brower.display(protein_name)
def folding_landscape_3d_plot2(protein_name,native_pdb_name,pdb_path,configure_file,RMSD_matrix_file_name):

    if pdb_path[-1]!='/':
        print(pdb_path)
        pdb_path=pdb_path+'/'
    x_coords = []
    y_coords = []
    z_coords = []
    pdb_label_list=[]
 
    print(protein_name+'_coords.cfg')
    coord_config_file = "/Users/congwanwang/Desktop/loudou.txt"

    need_load_flag=False
    pdb_coords_file_reader=open(coord_config_file)
    line=pdb_coords_file_reader.readline()
    pdb_name_list2=[]
    while line:
        #print('h1')
        line_list=line.split()
        if len(line_list)>=4:
#                 print('2')
            pdb_name_list2.append(line_list[0])
            x_coord = float(line_list[1])
            y_coord = float(line_list[2])
            #if x_coord>-0.248 or y_coord>0.15 or y_coord< -0.15 or x_coord<-0.348 :
             #   line=pdb_coords_file_reader.readline()
              #  continue

            x_coords.append(x_coord)
            y_coords.append(y_coord)
            z_coords.append(float(line_list[3]))
            line=pdb_coords_file_reader.readline()
            pdb_label_list.append(line_list[0]+' '+line_list[4])
        else:
            need_load_flag=True
            print('3')
            break

    brower=data_brower(x_coords,y_coords,z_coords,pdb_label_list)
    brower.display(protein_name)

def test():   
#     configure_file='/home/kdy/Seafile/workspace/EnergyLandscape3DTest/src/FileBatchDeal/test_score_configure'
    configure_file='/Users/congwanwang/Desktop/loudou.txt'
#     configure_file='/home/kdy/Seafile/workspace/EnergyLandscape3DTest/src/experiments_data/rosetta_decoys_62proteins/rosetta1/1a32/1a32.score'
    pdb_name='/Users/congwanwang/Desktop'
    native_pdb_name='1ctf.pdb'
    pdb_path='/Users/congwanwang/Desktop'
#     pdb_path='/home/kdy/Seafile/workspace/EnergyLandscape3DTest/src/pdbtest/'
#     pdb_path='/home/kdy/Seafile/workspace/EnergyLandscape3DTest/src/experiments_data/rosetta_decoys_62proteins/rosetta1/1a32/'
    RMSD_matrix_file_name='RMSD_matrix'
    folding_landscape_3d_plot2(pdb_name,native_pdb_name,pdb_path,configure_file,RMSD_matrix_file_name)

test()

       
    
