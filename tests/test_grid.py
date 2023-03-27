import pytest
from pykasso import Grid


@pytest.fixture
def empty_grid():
    # TODO
    return Grid(0, 0, 0, 0, 0, 0, 0, 0, 0)


@pytest.fixture
def example_grid():
    # TODO
    return Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)


def test_grid_get_i(example_grid):
    assert example_grid.get_i(70) == 7
    
    
def test_grid_get_j(example_grid):
    assert example_grid.get_j(70) == 3
    
    
def test_grid_get_k(example_grid):
    assert example_grid.get_k(70) == 2


def test_grid_get_indices(example_grid):
    assert example_grid.get_indices(70, 70, 70) == (7, 3, 2)
    

def test_grid_get_x(example_grid):
    assert example_grid.get_x(5) == 50
    
    
def test_grid_get_y(example_grid):
    assert example_grid.get_y(5) == 100
    
    
def test_grid_get_z(example_grid):
    assert example_grid.get_z(5) == 150
    

def test_grid_get_coordinates(example_grid):
    assert example_grid.get_coordinates(5, 5, 5) == (50, 100, 150)
    
    
def test_grid_is_inbox(example_grid):
    assert example_grid.is_inbox(70, 70, 70)
    assert not example_grid.is_inbox(70, 70, 1000)
