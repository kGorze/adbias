from pathlib import Path
import math



class Draw():
    def __init__(self, draw_candidate_points, draw_rejected_points, draw_current_cube, draw_corrected_box, accepted_point_radius):
        self.draw_candidate_points = draw_candidate_points
        self.draw_rejected_points = draw_rejected_points
        self.draw_current_cube = draw_current_cube
        self.draw_corrected_box = draw_corrected_box
        self.accepted_point_radius = accepted_point_radius

def parse_autodock_mapfile(mapfile_path):
    spacing = None
    nelements = None
    center = None

    with open(mapfile_path, 'r') as file:
        for _ in range(6):
            line = file.readline()


class Bias():
    def __init__(self, mapfile_path, bias_list[bias_x,bias_y,bias_z], Vset, r, epsilon):
        self.map_path = Path(mapfile_path)
        self.bias_x = bias_list[0] if len(bias_list) > 0 else 0
        self.bias_y = bias_list[1] if len(bias_list) > 1 else 0
        self.bias_z = bias_list[2] if len(bias_list) > 2 else 0
        self.bias_center = (self.bias_x, self.bias_y, self.bias_z)
        if Vset >= 0:
            raise ValueError("Vset must be negative for attractive bias")
        self.Vset = Vset
        if(r <= 0):
            raise ValueError("r must be positive for attractive bias")
        self.r = r
        if(epsilon <= 0):
            raise ValueError("epsilon must be greater than zero for a finite numerical cutoff")
        self.epsilon = epsilon

    pass

# zwrocenie liczby miedzy wartościami, ukrucenie jej
def clamp_int(value, min_value, max_value):
    return max(min_value, min(value, max_value))

#indeksy w gridzie są nienegatywne. floor(u+0.5) daje deterministyczny indeks, jeżeli taka jest zasada. Zasada floor(u + 0.5) zaokrągla nieujemną liczbę do najbliższej liczby całkowitej, a dokładne połówki, np. 6.5, zawsze zaokrągla w górę.
def nearest_grid_index(coord, minimum, spacing, n_points):
    u = (coord - minimum) / spacing   
    idx = int(math.floor(u + 0.5))
    return clamp_int(idx, 0, n_points - 1)

# dystans euklidesowski 3D między dwoma punktami (x,y,z)
def distance3(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    return distance

def add3(a, b):
    x = a[0] + b[0]
    y = a[1] + b[1]
    z = a[2] + b[2]
    return (x, y, z)

def fracion_color():
    #make a gradient that just spit up the numbers to the tlc from red to yellow with less opacity. output as (fraction, color, opacity)
    #red = (1, (255,255,255), 0)

    #orange = (0.5, (255,165,0), 0.)

    #yellow = (0.0, (255, 255, 0), 0.0)
