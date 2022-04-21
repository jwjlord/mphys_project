# visualisations

# Here are some functions which use tk to make a pyrochlore lattice animation.
# they could easily be used to make animations of other lattices.
# these functions need some kind of tests writing for them

import tkinter as tk
from .exchange import *
from itertools import *

class DrawPyroLattice(tk.Canvas):
    """Facilities for drawing 3D diagrams of the pyrochlore lattice."""
    
    def __init__(self, spin_lattice, site_lattice, width=1000, height=1000, scale=500):
        
        import numpy as np
        
        self.frame=tk.Frame()
        tk.Canvas.__init__(self, master=self.frame, width=width, height=height)
        
        self.spin_lattice = spin_lattice
        self.site_lattice = site_lattice
        
        self.theta = 0
        self.phi = 0
        self.scale = scale
        
    def draw(self):
        """Displays the tkinter canvas.
        :return: None"""
        
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)
        
        self.bind("<Left>", self._left)
        self.bind("<Right>", self._right)
        self.bind("<Up>", self._up)
        self.bind("<Down>", self._down)
        self.bind("<Button-1>", self._zoom_in)
        self.bind("<Button-2>", self._zoom_out)

        self.grid(column=0, row=0)
        self.frame.grid(column=0, row=0)
        self.mainloop()
        
    def _up(self, event):
        import numpy as np
        self.theta += 0.01 * np.pi
        self.delete('all')
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)

    def _down(self, event):
        import numpy as np
        self.theta -= 0.01 * np.pi
        self.delete('all')
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)

    def _right(self, event):
        import numpy as np
        self.phi += 0.01 * np.pi
        self.delete('all')
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)

    def _left(self, event):
        import numpy as np
        self.phi -= 0.01 * np.pi
        self.delete('all')
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)

    def _zoom_in(self, event):
        self.scale = self.scale * 1.1
        self.delete('all')
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)

    def _zoom_out(self, event):
        self.scale = self.scale * (10/11)
        self.delete('all')
        self.draw_pyro_lattice(self.spin_lattice, self.site_lattice, self.theta, self.phi, self.scale)
        
    def draw3d(self, start_point, end_point, theta, phi, scale, offset):
        """Draw a 3D line on the 2D canvas
        :start_point: np.ndarray(float), starting point in 3D
        :end_point: np.ndarray(float), end point in 3D
        :theta: float, observer polar angle
        :phi: float, observer azimuthal angle
        :scale: float, scale factor
        :offset: float, image top left corner offset amount from canvas top left corner (equal in x and y directions)
        :return: None"""
        
        import numpy as np
        
        start_end = [start_point - offset, end_point - offset]
        
        start_end_trans = [
            scale * (rotate(point, theta, phi) + offset)
            for point in start_end
        ]
        
        self.create_line(start_end_trans[0][0], start_end_trans[0][1], start_end_trans[1][0], start_end_trans[1][1])
        
    def draw3dpoly(self, theta, phi, scale, offset, color, c1, c2, c3, *args):
        """Draw a 3D polygon on the 2D canvas
        :theta: float, observer polar angle
        :phi: float, observer azimuthal angle
        :scale: float, scale factor
        :offset: float, image first corner offset amount from canvas top left corner (equal in x and y directions)
        :color: str, tkinter colour to fill the face
        :c1: np.ndarray(float), corner 1 coordinates
        :c2: np.ndarray(float), corner 1 coordinates
        :c3: np.ndarray(float), corner 1 coordinates
        :args: np.ndarray(float), additional corner coordinates
        :return: None"""
        
        from itertools import cycle
        import numpy as np
        
        corners = [
            scale * (
                rotate(corner - offset, theta, phi) + offset
            )
            for corner in [c1, c2, c3] + list(args)
        ]
        
        corners_2d = [
            corners[i][j]
            for i, j in chain.from_iterable(
                zip(
                    repeat(i, 2),
                    range(0,2)
                ) for i in range(0, len(corners))
            )
        ]


        self.create_polygon(corners_2d, fill=color)
        
    def arrow(self, start_point, end_point, theta, phi, scale, offset, color):
        """Draw a 3D arrow on the 2D canvas. Arrow points towards end_point.
        :start_point: np.ndarray(float), starting point in 3D
        :end_point: np.ndarray(float), end point in 3D
        :theta: float, observer polar angle
        :phi: float, observer azimuthal angle
        :scale: float, scale factor
        :offset: float, image top left corner offset amount from canvas top left corner (equal in x and y directions)
        :color: str, arrow tip colour
        :return: None"""
        
        import numpy as np
        
        direction = end_point - start_point # 3-vector direction
        
        # arrow shaft
        
        self.draw3d(
            start_point + 0.2 * direction, # shorten arrow slightly
            end_point - 0.2 * direction,
            theta,
            phi,
            scale,
            offset
        )
        
        # Now locate the arrow tips using perpendiculars to the shaft
        
        arrow_len = np.linalg.norm(direction)
        direction_norm = direction / arrow_len
        
        # first: an arbitrary perpendicular
        
        p1 = get_perp(direction_norm)
        
        # then the mutual perpendicular
        
        p2 = np.cross(direction_norm, p1)
        
        # draw two triangles for the arrow tip
        
        for p in [p1, p2]:
        
            self.draw3dpoly(
                theta,
                phi,
                scale,
                offset,
                color,
                end_point - 0.2 * direction,
                end_point - 0.3 * direction + 0.2 * arrow_len * p,
                end_point - 0.3 * direction - 0.2 * arrow_len * p,
            )
            
    def draw_text(self, value, location, theta, phi, scale, offset):
        """Draw text on the canvas.
        :value: str, the text to be drawn
        :location: coordinates to draw text at
        :theta: float, observer polar angle
        :phi: float, observer azimuthal angle
        :scale: float, scale factor
        :offset: float, top left corner offset amount from canvas top left corner (equal in x and y directions)
        :return: None"""
        
        location_trans = scale * (
            rotate(location - offset, theta, phi) + offset
        )
        
        self.create_text(location_trans[0], location_trans[1], text=f'{value}')
        
        
    def draw_cube(self, tl, scale, side, offset, theta, phi):
        """Draws a cube
        :tl: np.ndarray(float), coordinates of front top left corner
        :scale: float, scale factor
        :side: float, side length
        :offset: float, top left corner offset amount from canvas top left corner (equal in x and y directions)
        :theta: float, observer polar angle
        :phi: float, observer azimuthal angle
        :return: None"""
        
        import numpy as np
        
        # Get the eight cube vertices
        # We make 3-tuples with all possible combinations of 0 and side, then collect all the unique permutations of these. This gives the vertices.
        
        cube_vertices = {
            *chain.from_iterable(
                permutations(i, 3)
                for i in combinations_with_replacement([0, side], 3)
            )
        }
        
        # convert to arrays. can't be done in previous step as numpy array is not hashable
        
        cube_vertices_arrays = [np.array(j) + tl for j in cube_vertices]
        
        for element in permutations(cube_vertices_arrays, 2):
            
            if np.linalg.norm(element[1] - element[0]) == side:
                
                self.draw3d(element[0], element[1], theta, phi, scale, offset)
            
            else:
                pass

        
    def draw_pyro_lattice(self, spin_lattice, site_lattice, theta, phi, scale):
        
        import numpy as np
        
        diam_basis = np.array([[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.25,0.25,0.25],[0.25,0.75,0.75],[0.75,0.25,0.75],[0.75,0.75,0.25]])
        
        ice_basis = np.array([[0,0,0],[0,0,0],[0,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0.5,0],[0.25,0.25,0.25],[0.25,0.25,0.25],[0.25, 0.75, 0.75],[0.25, 0.75, 0.75],[0.75, 0.25, 0.75],[0.75, 0.25, 0.75],[0.75, 0.75, 0.25],[0.75, 0.75, 0.25]])
        
        side = site_lattice.shape[0]
        a = 0.5
        offset = 0.5 * side + a
        
        # draw two faces with colour

        self.draw3dpoly(
            theta,
            phi,
            scale,
            offset,
            'IndianRed1',
            np.array([a, a, a]),
            np.array([a, a+side, a]),
            np.array([a+side, a+side, a]),
            np.array([a, a+side, a])
        )
        
        self.draw3dpoly(
            theta,
            phi,
            scale,
            offset,
            'CadetBlue1',
            np.array([a, a+side, a]),
            np.array([a, a+side, a+side]),
            np.array([a, a, a]),
            np.array([a, a, a+side])
        )
        
        # Make the cube outline
        
        self.draw_cube(
            np.array([a, a, a]),
            scale,
            side,
            offset,
            theta,
            phi
        )

        # Iterate over site lattice and place site values on the sites
        
        site_iter = np.nditer(site_lattice, op_flags = ['readwrite'], flags = ['multi_index'])
        
        for site in site_iter:
            
            # Get site location in real space
            
            site_index = np.array(site_iter.multi_index)
            site_loc = site_index[0:3] + diam_basis[(site_index[3], ...)] + a
            
            # periodic boundary conditions
            # if any of the sites are at the edge of the lattice:
            # fill in the corresponding sites on the opposite side
            
            if np.isclose(a, site_loc).any():
                
                # these give all possible displacements:
                # unique permutations of (0, side) with length 3
                
                displacements = {
                    *chain.from_iterable(
                        permutations(i, 3)
                        for i in combinations_with_replacement([0, side], 3)
                    )
                }
                
                # filter possible displacements to affect only the relevant dimensions
                
                for element in displacements:
                
                    site_loc_2 = site_loc + (
                        np.array(element) * np.isclose(np.array([a, a, a]), site_loc)
                    )
                
                    self.draw_text(
                        site - 2,
                        site_loc_2,
                        theta,
                        phi,
                        scale,
                        side
                    )
            
            # Mark number of monopole excitations on that location with text
            
            self.draw_text(
                site - 2,
                site_loc,
                theta,
                phi,
                scale,
                side
            )

        # Iterate over spin lattice and place spins as arrows between sites
        
        spin_iter = np.nditer(spin_lattice, op_flags = ['readwrite'], flags = ['multi_index'])
        
        for spin in spin_iter:
            
            spin_index = np.array(spin_iter.multi_index)
            
            # make an arrow
            
            self.place_arrow(spin, spin_index, theta, phi, scale, offset, a, side)

            # Periodic boundary conditions:
            
            if np.isclose(0, spin_index).any():
                
                # these give all possible displacements:
                # unique permutations of (0, side) with length 3
                
                displacements = {
                    *chain.from_iterable(
                        permutations(i, 3)
                        for i in combinations_with_replacement([0, side], 3)
                    )
                }
                
                # filter possible displacements to affect only the relevant dimensions
                # the place_arrow function will remove any arrows which sit outside the lattice
                
                for element in displacements:
                
                    spin_index_2 = spin_index + np.append(
                        np.array(element) * np.isclose(np.array([0, 0, 0]), spin_index[0:3]),
                        0
                    )
                
                    self.place_arrow(spin, spin_index_2, theta, phi, scale, offset, a, side)



    def place_arrow(self, spin, spin_loc, theta, phi, scale, offset, a, side):
        """Make an arrow representing a spin in the pyrochlore lattice.
        :spin: int, value of the spin (1 or -1)
        :spin_loc: index location of the spin
        :theta: observer polar angle
        :phi: observer azimuthal angle
        :scale: float, scale factor
        :offset: float, offset from top left corner
        :a: float, 
        :side: float, lattice side length. Required to handle periodic boundary conditions"""
        
        # spin to site key, as described in exchange energy functions
        
        import numpy as np
        
        diam_basis = np.array([[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.25,0.25,0.25],[0.25,0.75,0.75],[0.75,0.25,0.75],[0.75,0.75,0.25]])
        
        key = np.array([
            [[0,0,0,0],[0,-1,-1,5]],
            [[0,0,0,0],[-1,0,-1,6]],
            [[0,0,0,1],[0,0,0,5]],
            [[0,0,0,1],[-1,0,0,7]],
            [[0,0,0,2],[0,0,0,4]],
            [[0,0,0,2],[0,0,0,6]],
            [[0,0,0,3],[0,0,0,4]],
            [[0,0,0,3],[0,0,0,7]],
            [[0,0,0,4],[0,0,0,0]],
            [[0,0,0,4],[0,0,0,1]],
            [[0,0,0,5],[0,1,0,2]],
            [[0,0,0,5],[0,0,1,3]],
            [[0,0,0,6],[1,0,0,1]],
            [[0,0,0,6],[0,0,1,3]],
            [[0,0,0,7],[1,1,0,0]],
            [[0,0,0,7],[0,1,0,2]],
        ])
        
        site_loc = site_to_lattice(spin_loc, key)
        
        # start and end points
        
        points = [
            sl[0:3] + diam_basis[(sl[3], ...)] + a
            for sl in site_loc
        ]
        
        # make sure spin is inside the lattice
        
        condition = (
            (
                np.isclose(points[0], np.array([a + side, a + side, a + side])) | (points[0] < np.array([a + side, a + side, a + side]))
            ).all() and (
                np.isclose(points[0], np.array([a, a, a])) | (points[0] > np.array([a, a, a]))
            ).all() and (
                np.isclose(points[1], np.array([a + side, a + side, a + side])) | (points[1] < np.array([a + side, a + side, a + side]))
            ).all() and (
                np.isclose(points[1], np.array([a, a, a])) | (points[1] > np.array([a, a, a]))
            ).all()
        )
        
        if spin == 1 and condition:
            self.arrow(points[0], points[1], theta, phi, scale, offset, 'Black')
        
        elif spin == -1 and condition:
            self.arrow(points[1], points[0], theta, phi, scale, offset, 'Black')
            
        elif spin not in [1, -1]:
            raise ValueError('Lattice contains invalid spin value')
            
        else:
            pass
        
        
def rotate(vector, theta, phi):
    """Rotate a 3-vector in 3D polar coordinates
    :vector: np.ndarray(float), 3-vector to rotate
    :theta: polar rotation angle (radians)
    :phi: azimuthal rotation angle (radians)"""
    
    import numpy as np
    
    rot = np.array(
        [[np.cos(phi), 0, np.sin(phi)],
        [np.sin(theta) * np.sin(phi), np.cos(theta), -np.sin(theta) * np.cos(phi)],
        [-np.cos(theta) * np.sin(phi), np.sin(theta), np.cos(theta) * np.cos(phi)]]
    )
    
    return np.matmul(vector, rot)

def get_perp(vector):
    """Get an arbitrary perpendicular to a normal 3-vector
    :vector: np.ndarray(float), normalised 3-vector
    :return: np.ndarray(float), a normal perpendicular to the 3-vector"""
    
    import numpy as np
    
    perp = np.array(
        [
            0,
            - vector[2],
            vector[1]
        ]
    ) / np.sqrt(1 - vector[0]**2)
    
    return perp