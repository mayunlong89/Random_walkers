# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 14:26:50 2020

@author: Yunlong Ma
@E-mail: glb-biotech@zju.edu.cn

"""

#randomwalk_visualization.py 

import matplotlib.pyplot as plt
import random

from random_walk import RandomWalk

while True:
    #establish an example for RandWalk
    rw = RandomWalk(50000)
    rw.fill_walk()
    
    #set the size of drawing window
    plt.figure(figsize=(10,6))
    
    #mapping
    point_numbers = list(range(rw.num_points))
    plt.scatter(rw.x_values,rw.y_values, edgecolors="none",s=1,c=point_numbers,cmap=plt.cm.Blues)
    
    #highlight starting and ending points
    plt.scatter(0,0,c="green",edgecolors="none",s=5)
    plt.scatter(rw.x_values[-1],rw.y_values[-1],c='red',edgecolors='none',s=5)
    
    #labelling
    plt.xlabel('Axis X', fontsize=14)
    plt.ylabel('Axis Y', fontsize=14)
    plt.title("Random Walk Visualization",fontsize=24)
    
    #Hidden axes
    plt.axes().get_xaxis().set_visible(False)
    plt.axes().get_yaxis().set_visible(False)
    
    #visualization
    plt.show()  #this code can not be used with save figure in the same time, otherwise, the figure saved with blank
    
    #save figure
    #filename = 'rw_' + str(random.randint(000000,999999)) + '.png'
    #plt.savefig('images/'+filename, bbox_inches='tight')
    
    #ask contiune or not???
    flag = input('Make another walk? (y/n):')
    if flag.lower() not in ['y','yes']:
        break
    
    
    
