#2020-08-18
#@Author: Yunlong Ma
#@E-mail: glb-biotech@zju.edu.cn

#random_walk.py

from random import choice

class RandomWalk():
    '''
    Generate a class of Random walk for data
    '''
    def __init__(self, num_points=5000):
        '''
        Initialization parameters

        Parameters
        ----------
        num_points : TYPE, optional
            DESCRIPTION. The default is 5000.

        Returns
        -------
        None.

        '''
        self.num_points = num_points
        
        #all random walks start by (0,0)
        self.x_values = [0]
        self.y_values = [0]
    
    def fill_walk(self):
        '''
        Calculate all points contained by random walks
        '''
        
        #with random walk goes, the given length of list should be setted
        while len(self.x_values) < self.num_points:
            
            #deside the direction and distance for random walk
            x_direction = choice([1,-1])
            x_distance = choice([0,1,2,3,4])
            x_step = x_direction*x_distance
            
            y_direction = choice([1,-1])
            y_distance = choice([0,1,2,3,4])
            y_step = y_direction*y_distance
            
            #reject march on the spot
            if x_step == 0 and y_step == 0:
                continue
            
            #calculate the x and y values of next point
            next_x = self.x_values[-1] + x_step
            next_y = self.y_values[-1] + y_step
            
            self.x_values.append(next_x)
            self.y_values.append(next_y)
    
        

