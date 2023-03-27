"""
TODO
"""

### Internal dependencies
import os
import pickle

### External dependencies
import yaml
import numpy as np
import pandas as pd


class ProjectReader():
    """
    TODO
    """
    
    def __init__(self, project_directory: str, *args, **kwargs) -> None:
        # Reads the project_state.yaml file
        self.project_state = self.open_project(project_directory)
        
    def open_project(self, path: str) -> dict:
        path = path.strip("/")
        path = path + "/outputs/project_state.yaml"
        with open(path, "r") as f:
            return yaml.safe_load(f)
    
    def _get_grid_origin(self) -> tuple:
        x0 = self.project_state['grid']['x0']
        y0 = self.project_state['grid']['y0']
        z0 = self.project_state['grid']['z0']
        return (x0, y0, z0)
    
    def _get_grid_dimensions(self) -> tuple:
        nx = self.project_state['grid']['nx']
        ny = self.project_state['grid']['ny']
        nz = self.project_state['grid']['nz']
        return (nx, ny, nz)
    
    def _get_grid_spacing(self) -> tuple:
        dx = self.project_state['grid']['dx']
        dy = self.project_state['grid']['dy']
        dz = self.project_state['grid']['dz']
        return (dx, dy, dz)

    def _read_pickle(self, path: str):
        with open(path, 'rb') as handle:
            return pickle.load(handle)


class GSLIB():
    """
    TODO
    """
    
    def __init__(self, path):
        if not os.path.exists(path):
            msg = "The path is not valid. '{}' does not exist.".format(path)
            raise FileNotFoundError(msg)
        else:
            self.location = path

        with open(path) as f:
            lines = f.readlines()
            
        self.comment = lines[0].strip()
        
        try:
            self.vars_nbr = int(lines[1])
        except Exception as err:
            raise err
            
        self.vars_names = []
        for i in range(self.vars_nbr):
            self.vars_names.append(lines[i + 2].strip())
        
        try:
            skip = 2 + self.vars_nbr
            data = np.genfromtxt(path, skip_header=skip)
            self.data = pd.DataFrame()
            if self.vars_nbr == 1:
                self.data[self.vars_names[0]] = data
            else:
                for i, name in zip(range(self.vars_nbr), self.vars_names):
                    self.data[name] = data[:, i]
        except Exception as err:
            raise err
        
    def __str__(self):
        txt = "GSLIB File Reader \n"
        txt += "File location         : {}\n".format(self.location)
        txt += "Comment line          : {}\n".format(self.comment)
        txt += "Number of variable(s) : {}\n".format(self.vars_nbr)
        txt += "Names of variable(s)  : {}\n".format(self.vars_names)
        txt += "Length of data        : {}".format(len(self.data))
        return txt
     
    def export_column(self, column, dim):
        
        data = self.data.loc[:, column].values
        
        if len(dim) == 2:
            nx, ny = dim
            data = np.reshape(data, (nx, ny), order='F')
        
        elif len(dim) == 3:
            nx, ny, nz = dim
            data = np.reshape(data, (nx, ny, nz), order='F')
            
        return data
