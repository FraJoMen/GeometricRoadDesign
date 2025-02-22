# -*- coding: utf-8 -*-
"""


@author: Francisco 
"""


#%%
import numpy as np
import plotly.graph_objects as go
#import plotly.io as pio

#%%
class CircularTransition:
    """
    Class to handle the circular transition between two straight road sections.
    """
    def __init__(self, P0, P1, P2, R):
        self.P0 = np.array(P0)
        self.P1 = np.array(P1)
        self.P2 = np.array(P2)
        self.R = R
        
        # Compute unit direction vectors
        self.v1 = (self.P0 - self.P1) / np.linalg.norm(self.P1 - self.P0)
        self.v2 = (self.P2 - self.P1) / np.linalg.norm(self.P2 - self.P1)
        
        self.alpha_deg = self.compute_intersection_angle()
        self.theta_deg = 180 - self.alpha_deg  # θ = π - α in degrees
        
        # Compute characteristic lengths
        self.T = self.R * np.tan(np.radians(self.theta_deg) / 2)
        self.B = self.R * np.tan(np.radians(self.alpha_deg) / 2)
        
        self.O = self.compute_circle_center()
        self.T1, self.T2 = self.compute_tangent_points()
        
        # Print parameters immediately after initialization
        self.print_parameters()
          
    def print_parameters(self):
        """
        Prints the geometric parameters of the circular transition.
        """
        print(f"P0: {self.P0}, P1: {self.P1}, P2: {self.P2}")
        print(f"Radius (R): {self.R}")
        print(f"Intersection Angle α: {self.alpha_deg:.2f} degrees")
        print(f"Central Angle θ: {self.theta_deg:.2f} degrees")
        print(f"Tangent Length T: {self.T:.4f}")
        print(f"Bisector Length B: {self.B:.4f}")
        print(f"Circle Center (O): {self.O}")
        print(f"Tangent Points: T1 = {self.T1}, T2 = {self.T2}")
      
    def compute_intersection_angle(self):
        """
        Computes the intersection angle between two straight road segments.
        """
        dot_product = np.dot(self.v1, self.v2)
        alpha_rad = np.arccos(dot_product)
        return np.degrees(alpha_rad)
    
    def compute_circle_center(self):
        """
        Computes the center of the circular transition using perpendicular bisector method.
        """
        bisector = (self.v1 + self.v2)
        bisector /= np.linalg.norm(bisector)
        
        O = self.P1 + bisector * (self.B + self.R)  # Corrected formula using bisector length B and radius R
        return O
    
    def compute_tangent_points(self):
        """
        Computes the tangent points T1 and T2 where the straight segments meet the circle.
        """
        T1 = self.P1 + self.v1 * self.T  # Corrected T1 using tangent length T
        T2 = self.P1 + self.v2 * self.T  # Corrected T2 using tangent length T
        
        return T1, T2
    
    def plot(self):
        """
        Plots the transition with relevant geometric elements.
        """
        x_coords = [self.P0[0], self.P1[0], self.P2[0], self.O[0], self.T1[0], self.T2[0]]
        y_coords = [self.P0[1], self.P1[1], self.P2[1], self.O[1], self.T1[1], self.T2[1]]
        
        fig = go.Figure()
        
        # Adding points
        fig.add_trace(go.Scatter(x=x_coords[:-3], y=y_coords[:-3], mode='markers+text',
                                 marker=dict(size=10, color='red'),
                                 text=['P0', 'P1', 'P2'], textposition="top center",
                                 name='Points'))
        
        # Adding circle center
        fig.add_trace(go.Scatter(x=[self.O[0]], y=[self.O[1]], mode='markers+text',
                                 marker=dict(size=10, color='blue'),
                                 text=['O'], textposition="top center",
                                 name='Circle Center'))
        
        # Adding tangent points
        fig.add_trace(go.Scatter(x=[self.T1[0], self.T2[0]], y=[self.T1[1], self.T2[1]],
                                 mode='markers+text', marker=dict(size=10, color='green'),
                                 text=['T1', 'T2'], textposition="top center",
                                 name='Tangent Points'))
        
        # Adding segment lines
        fig.add_trace(go.Scatter(x=x_coords[:-3], y=y_coords[:-3], mode='lines',
                                 line=dict(dash='dash', color='black'),
                                 name='Segments'))
        
        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        
        fig.update_layout(
            title="Interactive Transition Plot", 
            xaxis_title="X", yaxis_title="Y", 
            template="plotly_white",
            legend=dict(orientation="h", yanchor="bottom", y=-0.4, xanchor="center", x=0.5)
        )
        
        fig.show()

# Example usage (Uncomment to use in a script)
# transition = CircularTransition((0, 0), (2, 2), (4, 0), 5)
# transition.plot()
