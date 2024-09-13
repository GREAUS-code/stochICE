
"""
Class to calculate hydraulic parameters for a normal depth steady-state flow 
in an irregular cross-section. This class helps construct a stage-discharge
curve adapted in the presence of an ice cover at the downstream end of the model. 
This code is only used however to obtain the wetted area (A) and wetted perimeter (P) 
for a given water level. The open-channel discharge this code produces is discarded. 
The discharge for an ice covered channel is obtained in stochRIVICE.py by solving 
the Manning's equation for Q, knowing A and P+ice cover width, obtained from this
code. 
"""

import math

class IrregularSection:
    
    # This class is credited to https://github.com/alexiusacademia/ChannelFlowLib

    
    def __init__(self, points):
        """
        Constructor and initializations
        :param points:
        :return:
        """
        self.points = points
        # Initializations
        self.roughness = 0.0            # Average roughness of the cross section
        self.bed_slope = 0.0            # River bed average slope
        self.water_elevation = 0.0      # Water surface elevation
        self.wetted_area = 0.0          # Wetted area
        self.wetted_perimeter = 0.0     # Wetted perimeter
        self.hydraulic_radius = 0.0     # Hydraulic radius
        self.velocity = 0.0             # Average velocity
        self.discharge = 0.0            # Discharge
        self.max_water_elevation = 0.0
        self.min_water_elevation = 0.0
        self.froude_number = 0.0
        self.state=True

    # ---------
    # Setters
    # ---------
    def set_average_rougness(self, roughness_coefficient):
        """
        Set the manning's roughness coefficient, n
        :param roughness_coefficient:
        :return:
        """
        self.roughness = roughness_coefficient

    def set_bed_slope(self, bed_slope):
        """
        Set the average bed slope
        :param bed_slope:
        :return:
        """
        self.bed_slope = bed_slope

    def set_water_elevation(self, water_elevation):
        """
        Sets the water surface elevation
        :param water_elevation:
        :return:
        """
        self.water_elevation = water_elevation

    # ----------
    # Methods
    # ----------
    def analyze(self):
        # Validate inputs
        if self.bed_slope == 0:
            raise Exception
        if self.roughness == 0:
            raise Exception

        # Count the points
        num_points = len(self.points)

        # Get the lower bank elevation
        if self.points[0][1] > self.points[num_points-1][1]:
            max_ws = self.points[num_points-1][1]
        else:
            max_ws = self.points[0][1]

        self.max_water_elevation = max_ws

        if self.water_elevation > max_ws:
            print('Water will overflow the bank!')
            # raise Exception
            self.state=False
            return

        # Get the lowest possible water elevation
        if self.water_elevation < self.get_lowest_elev(self.points):
            print('Water surface is below the lowest point of the channel.')
            self.state=False
            return

        # Number of intersections
        left = 0
        right = 0

        new_points = []     # New points, removing points out of range of the intersection

        for index in range(len(self.points)):
            x, y = self.points[index]

            # Look for the first point of intersection
            if left == 0:
                if y < self.water_elevation:
                    left += 1
                    x1 = self.points[index-1][0]
                    y1 = self.points[index-1][1]
                    x2 = self.points[index][0]
                    y2 = self.points[index][1]
                    x3 = (self.water_elevation - y1) * (x2 - x1) / (y2 - y1) + x1
                    new_points.append((x3, self.water_elevation))  # Add this point as the first new point from the left

            if right == 0:
                if left == 1:
                    if y > self.water_elevation:
                        right += 1
                        x1 = self.points[index-1][0]
                        y1 = self.points[index-1][1]
                        x2 = self.points[index][0]
                        y2 = self.points[index][1]
                        x3 = (self.water_elevation - y1) * (x2 - x1) / (y2 - y1) + x1
                        new_points.append((x3, self.water_elevation))

            if left == 1:
                if right == 0:
                    new_points.append(self.points[index])

        left_top_bank_point = new_points[0]
        right_top_bank_point = new_points[len(new_points) - 1]

        # Hydraulic elements
        self.wetted_area = self.polygon_area(new_points)
        self.wetted_perimeter = self.get_perimeter(new_points)
        self.hydraulic_radius = self.wetted_area / self.wetted_perimeter
        self.velocity = (1 / self.roughness) * self.hydraulic_radius**(2/3) * self.bed_slope**0.5
        self.discharge = self.velocity * self.wetted_area

        self.top_width = right_top_bank_point[0] - left_top_bank_point[0]
        hydraulic_depth = self.polygon_area(new_points) / self.top_width
        self.froude_number = self.velocity / math.sqrt(9.81 * hydraulic_depth)
        self.discharge_intensity = self.discharge / self.top_width

        self.state=True

    def polygon_area(self, vertices):
        """
        Implementation of Shoelace Formula in finding the area of a closed
        polygon bounded by vertices
        :param vertices:
        :return:
        """
        n = len(vertices) # of corners
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += vertices[i][0] * vertices[j][1]
            area -= vertices[j][0] * vertices[i][1]
        area = abs(area) / 2.0
        return area

    def get_perimeter(self, points):
        """
        Get the total distance covered by multiple points
        :param points:
        :return:
        """
        p = 0.0         # perimeter
        n = len(points)
        for i in range(n-1):
            p1 = points[i]
            p2 = points[i+1]
            p += self.point_distance([p1, p2])

        return p

    def point_distance(self, points):
        """
        Get the distance between two points
        :param points:
        :return:
        """
        p1 = points[0]
        p2 = points[1]
        x1, y1 = p1
        x2, y2 = p2

        dist = math.sqrt((y2-y1)**2 + (x2-x1)**2)

        return dist

    def get_lowest_elev(self, points):
        """
        Get the lowest point from the vertices
        :param points:
        :return: lowest
        """
        elevs = []                  # List of elevations (ordinates)
        lowest = 0                  # Initial value of lowest
        for point in points:
            elevs.append(point[1])  # Iterate through the points and collect the ordinates
        for i in range(len(elevs)): # Find the lowest in the list of ordinates
            if i == len(elevs):
                break
            if elevs[i] < lowest:
                lowest = elevs[i]
            else:
                lowest = lowest
        return lowest