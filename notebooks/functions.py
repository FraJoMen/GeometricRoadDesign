# -*- coding: utf-8 -*-
"""
GeometricRoadDesign - Road Design Functions
===========================================

This module is part of the GeometricRoadDesign project and provides tools 
for analyzing and designing road alignments, including:

1. **Planimetric Curve Transitions**
   - Circular transitions between straight road sections.
   - Future development: Clothoid-based transitions.

2. **Altimetric Curve Transitions (Vertical Alignments)**
   - Implementation planned.

3. **Complete Road Design Verification**
   - Future tools for checking geometric and regulatory compliance.

Current Implementation:
-----------------------
- `CircularTransition`: Computes the geometric parameters of a circular transition 
  (radii, intersection angles, tangent lengths, bisectors) and generates a visualization.

This module is designed for scalability, allowing future expansion with additional 
road design tools.

Author: Francisco
"""

import numpy as np
import plotly.graph_objects as go

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
        self.B = self.R / np.cos(np.radians(self.theta_deg) / 2)
        
        self.O = self.compute_circle_center()
        self.T1, self.T2 = self.compute_tangent_points()

        # Compute arc directly in global coordinates
        self.arc_x, self.arc_y = self.compute_arc_coordinates()

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
        
        O = self.P1 + bisector * self.B
        return O
    
    def compute_tangent_points(self):
        """
        Computes the tangent points T1 and T2 where the straight segments meet the circle.
        """
        T1 = self.P1 + self.v1 * self.T  
        T2 = self.P1 + self.v2 * self.T  
        return T1, T2

    def compute_arc_coordinates(self):
        """
        Computes and stores the coordinates of the circular arc in the global system.
        """
        # Compute start and end angles
        theta1 = np.arctan2(self.T1[1] - self.O[1], self.T1[0] - self.O[0])
        theta2 = np.arctan2(self.T2[1] - self.O[1], self.T2[0] - self.O[0])

        # Ensure correct arc direction
        theta_center = abs(theta2 - theta1)
        if theta_center > np.pi:
            if theta1 < theta2:
                theta1 += 2 * np.pi
            else:
                theta2 += 2 * np.pi

        theta_vals = np.linspace(theta1, theta2, 100)
        arc_x = self.O[0] + self.R * np.cos(theta_vals)
        arc_y = self.O[1] + self.R * np.sin(theta_vals)

        return arc_x, arc_y

    def plot(self):
        """
        Plots the transition with all elements: arc, segments, points, and vectors.
        """
        fig = go.Figure()

        # Add arc to plot (in blue)
        fig.add_trace(go.Scatter(x=self.arc_x, y=self.arc_y, mode='lines', name='Circular Arc', line=dict(color='blue', width=2)))

        # Add straight road segments P0P1 and P1P2 (black, dashed)
        fig.add_trace(go.Scatter(x=[self.P0[0], self.P1[0]], y=[self.P0[1], self.P1[1]], 
                                 mode='lines', line=dict(color='black', width=2, dash='dash'),
                                 name='P0P1'))
    
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P2[0]], y=[self.P1[1], self.P2[1]], 
                                 mode='lines', line=dict(color='black', width=2, dash='dash'),
                                 name='P1P2'))

        # Add reference segments P1O, T1O, T2O (black, dashed)
        fig.add_trace(go.Scatter(x=[self.O[0], self.P1[0]], y=[self.O[1], self.P1[1]], 
                                 mode='lines', line=dict(color='black', width=2, dash='dot'),
                                 name='P1O'))
        fig.add_trace(go.Scatter(x=[self.O[0], self.T1[0]], y=[self.O[1], self.T1[1]], 
                                 mode='lines', line=dict(color='black', width=2, dash='dot'),
                                 name='T1O'))
        fig.add_trace(go.Scatter(x=[self.O[0], self.T2[0]], y=[self.O[1], self.T2[1]], 
                                 mode='lines', line=dict(color='black', width=2, dash='dot'),
                                 name='T2O'))

        # Add reference points P0, P1, P2, O, T1, T2 (red markers)
        fig.add_trace(go.Scatter(
            x=[self.P0[0], self.P1[0], self.P2[0], self.O[0], self.T1[0], self.T2[0]],
            y=[self.P0[1], self.P1[1], self.P2[1], self.O[1], self.T1[1], self.T2[1]],
            mode='markers+text',
            marker=dict(size=10, color='red'),
            text=['P0', 'P1', 'P2', 'O', 'T1', 'T2'],
            textposition="top center",
            name='Reference Points'
        ))

        # Adjust plot size
        fig.update_layout(width=900, height=700, template="plotly_white")

        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)

        fig.update_layout(title="Circular Transition Plot", xaxis_title="X", yaxis_title="Y")

        fig.show()
