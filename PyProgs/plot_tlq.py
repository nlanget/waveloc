import os
import numpy as np
# import matplotlib.pyplot as plt
from visualization import setup_test_grid
    
from traits.api import *
from traitsui.api import *

from chaco.api import *

from enable.api import Component, ComponentEditor

class Viewer(HasTraits):
    pass

class Grid4D(HasTraits):
    grid4D = Instance(np.ndarray)
    grid_dict = Dict
    dims = List
    
    # probs = Enum
    
    low = Int(0)
    max = Int(1)
    current = Range('low','max')
    
    pd = ArrayPlotData()
    
   
    p_x0 = Plot(pd)
    p_x1 = Plot(pd)
    p_x2 = Plot(pd)
    p_x3 = Plot(pd)
    
    p_x0_x1 = Plot(pd)
    p_x0_x2 = Plot(pd)
    p_x0_x3 = Plot(pd)
    p_x1_x2 = Plot(pd)
    p_x1_x3 = Plot(pd)
    p_x2_x3 = Plot(pd)
    
    def __init__(self):
        self.grid4D,self.grid_dict, x_list = setup_test_grid()
        self.x0,self.x1,self.x2,self.x3 = x_list
        
        for key in self.grid_dict.keys():
            print key, self.grid_dict[key].shape
        
        self._load()
        
        self.p_x0.plot(('x0','p_x0'))
        self.p_x1.plot(('x1','p_x1'))
        self.p_x2.plot(('x2','p_x2'))
        self.p_x3.plot(('x3','p_x3'))
        
        self.p_x0_x1.contour_plot('p_x0_x1',type="poly",
                                            xbounds=(self.x0[0], self.x0[-1]),
                                            ybounds=(self.x1[0], self.x1[-1]))
        self.p_x0_x2.contour_plot('p_x0_x2',type="poly",
                                            xbounds=(self.x0[0], self.x0[-1]),
                                            ybounds=(self.x2[0], self.x2[-1]))
        self.p_x0_x3.contour_plot('p_x0_x3',type="poly",
                                            xbounds=(self.x0[0], self.x0[-1]),
                                            ybounds=(self.x3[0], self.x3[-1]))
        self.p_x1_x2.contour_plot('p_x1_x2',type="poly",
                                            xbounds=(self.x1[0], self.x1[-1]),
                                            ybounds=(self.x2[0], self.x2[-1]))
        self.p_x1_x3.contour_plot('p_x1_x3',type="poly",
                                            xbounds=(self.x1[0], self.x1[-1]),
                                            ybounds=(self.x3[0], self.x3[-1]))
        self.p_x2_x3.contour_plot('p_x2_x3',type="poly",
                                            xbounds=(self.x2[0], self.x2[-1]),
                                            ybounds=(self.x3[0], self.x3[-1]))
    
    def _grid4D_changed(self):
        self.dims = list(self.grid4D.shape)
        self.max = self.dims[-1]
    
    def _load(self):
        for var in ['x0','x1','x2','x3','x0_x1','x0_x2','x0_x3','x1_x2','x1_x3','x2_x3']:
            if len(var) > 2:
                var2 = var[:2]
            else:
                var2 = var
            self.pd.set_data(var,eval("self.%s"%var2))
            self.pd.set_data('p_%s'%var , self.grid_dict['prob_%s'% var])
    
    def _current_changed(self, old, new):
        self.pd.set_data('p_x0_x3i',np.ones(len(self.x0))*self.x3[new])
        self.p_x0_x3.plot(('x0','p_x0_x3i'))
        
        self.pd.set_data('p_x1_x3i',np.ones(len(self.x1))*self.x3[new])
        self.p_x1_x3.plot(('x1','p_x1_x3i'))
        
        self.pd.set_data('p_x2_x3i',np.ones(len(self.x2))*self.x3[new])
        self.p_x2_x3.plot(('x2','p_x2_x3i'))
        
        self.pd.set_data('p_x3i',np.ones(len(self.x2))*self.x3[new])
        self.p_x3.plot(('p_x3i', 'p_x3'))
        

    view = View("current",
            HGroup(
                Item("p_x0",editor=ComponentEditor()),
                Item("p_x1",editor=ComponentEditor()),
                Item("p_x2",editor=ComponentEditor()),
                Item("p_x3",editor=ComponentEditor()),
                ),    
            HGroup(
                Item("p_x0_x1",editor=ComponentEditor()),
                Item("p_x0_x2",editor=ComponentEditor()),
                Item("p_x0_x3",editor=ComponentEditor()),
                ),
            HGroup(
                Item("p_x1_x2",editor=ComponentEditor()),
                Item("p_x1_x3",editor=ComponentEditor()),
                Item("p_x2_x3",editor=ComponentEditor()),
                ),
                
                )
    
   




if __name__ == '__main__':
  g = Grid4D()
  g.configure_traits()
  

