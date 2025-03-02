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

class CircularTransition:
    """
    Class to handle the circular transition between two straight road sections.
    """
    def __init__(self, P0, P1, P2, R):
        self.P0 = np.array(P0)
        self.P1 = np.array(P1)
        self.P2 = np.array(P2)
        self.R = R
        
        # Ensure that the points are not collinear before proceeding
        if self.areCollinear(self.P0, self.P1, self.P2):
            raise ValueError("The three points are collinear. No circular transition is possible.")
        
        # Compute unit direction vectors
        self.v1 = (self.P0 - self.P1) / np.linalg.norm(self.P1 - self.P0)
        self.v2 = (self.P2 - self.P1) / np.linalg.norm(self.P2 - self.P1)
        
        self.alpha_deg = self.ComputeIntersectionAngle()
        self.theta_deg = 180 - self.alpha_deg  # θ = π - α in degrees
        
        # Compute characteristic lengths
        self.T = self.R * np.tan(np.radians(self.theta_deg) / 2)
        self.B = self.R / np.cos(np.radians(self.theta_deg) / 2)
        
        self.O = self.ComputeCircleCenter()
        self.T1, self.T2 = self.ComputeTangentPoints()
 
        # Compute arc directly in global coordinates
        self.arc_x, self.arc_y = self.ComputeArcCoordinates()

        # Print parameters immediately after initialization
        self.PrintParameters()

    @staticmethod
    def areCollinear(P0, P1, P2):
        """
        Determines if three points are collinear.
        """
        # Calculate the area of the triangle formed by the points
        area = 0.5 * np.linalg.norm(np.cross(P1 - P0, P2 - P0))
        return np.isclose(area, 0)

    @staticmethod
    def isInsideTriangle(P0, P1, P2, Pc):
        """
        Determines if the point Pc is inside the triangle formed by P0, P1, and P2.
        """
        def area(PA, PB, PC):
            return 0.5 * np.abs(PA[0] * (PB[1] - PC[1]) + PB[0] * (PC[1] - PA[1]) + PC[0] * (PA[1] - PB[1]))
        
        A = area(P0, P1, P2)
        A1 = area(Pc, P1, P2)
        A2 = area(P0, Pc, P2)
        A3 = area(P0, P1, Pc)
        
        return np.isclose(A, A1 + A2 + A3)
          
    def PrintParameters(self):
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
       
    def ComputeIntersectionAngle(self):
        """
        Computes the intersection angle between two straight road segments.
        """
        dot_product = np.dot(self.v1, self.v2)
        alpha_rad = np.arccos(dot_product)
        return np.degrees(alpha_rad)
    
    def ComputeCircleCenter(self):
        """
        Computes the center of the circular transition using perpendicular bisector method.
        """
        bisector = (self.v1 + self.v2)
        bisector /= np.linalg.norm(bisector)
        
        O = self.P1 + bisector * self.B
        return O
    
    def ComputeTangentPoints(self):
        """
        Computes the tangent points T1 and T2 where the straight segments meet the circle.
        """
        T1 = self.P1 + self.v1 * self.T  
        T2 = self.P1 + self.v2 * self.T  
        return T1, T2
 
    def ComputeArcCoordinates(self):
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
 
    def Plot(self, Pc=None, PA=None, PB=None):
        """
        Plots the transition with all elements: arc, segments, points, and vectors.
        If Pc is provided, it plots the control point as well.
        If PA and PB are provided, it plots the segment PA-PB as well.
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

        # Add control point Pc if provided
        if Pc is not None:
            fig.add_trace(go.Scatter(
                x=[Pc[0]],
                y=[Pc[1]],
                mode='markers+text',
                marker=dict(size=10, color='green'),
                text=['Pc'],
                textposition="top center",
                name='Control Point Pc'
            ))

        # Add segment PA-PB if provided
        if PA is not None and PB is not None:
            fig.add_trace(go.Scatter(
                x=[PA[0], PB[0]],
                y=[PA[1], PB[1]],
                mode='lines+markers+text',
                line=dict(color='orange', width=2),
                marker=dict(size=10, color='orange'),
                text=['PA', 'PB'],
                textposition="top center",
                name='Segment PA-PB'
            ))

        # Adjust plot size
        fig.update_layout(width=900, height=700, template="plotly_white")

        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)

        fig.update_layout(title="Circular Transition Plot", xaxis_title="X", yaxis_title="Y")

        fig.show()

    @staticmethod
    def FixPoint(P0, P1, P2, Pc):
        """
        Computes the unique radius R* such that the circular transition passes through the predefined point Pc.
        """
        P0 = np.array(P0)
        P1 = np.array(P1)
        P2 = np.array(P2)
        Pc = np.array(Pc)
        
        # Ensure that the points are not collinear before proceeding
        if CircularTransition.areCollinear(P0, P1, P2):
            raise ValueError("The three points are collinear. No circular transition is possible.")
        
        # Ensure that the point Pc is inside the triangle P0P1P2
        if not CircularTransition.isInsideTriangle(P0, P1, P2, Pc):
            raise ValueError("The point Pc is not inside the triangle formed by P0, P1, and P2.")
        
        # Compute unit direction vectors
        v1 = (P0 - P1) / np.linalg.norm(P1 - P0)
        v2 = (P2 - P1) / np.linalg.norm(P2 - P1)
        
        # Compute intersection angle α (in radians)
        alpha = np.arccos(np.dot(v1, v2))
        print(f"alpha (radians): {alpha}, alpha (degrees): {np.degrees(alpha)}")
        
        # Compute angle ω between vectors P1Pc and P1O
        P1Pc = Pc - P1
        P1O = (v1 + v2) / np.linalg.norm(v1 + v2)
        omega = np.arccos(np.dot(P1Pc / np.linalg.norm(P1Pc), P1O))
        print(f"omega (radians): {omega}, omega (degrees): {np.degrees(omega)}")
        
        # Compute angle φ at Pc in the triangle O P_c P_1
        phi = np.arcsin(np.sin(omega) / np.sin(alpha / 2))
        print(f"phi (radians): {phi}, phi (degrees): {np.degrees(phi)}")
        
        # Ensure φ is greater than 90 degrees
        if phi < np.pi / 2:
            phi = np.pi - phi
        
        # Compute angle γ
        gamma = np.pi - phi - omega
        print(f"gamma (radians): {gamma}, gamma (degrees): {np.degrees(gamma)}")
        
        # Compute the unique radius R*
        R_star = np.linalg.norm(Pc - P1) * np.sin(omega) / np.sin(gamma)
        print(f"R_star: {R_star}")
        
        # Return a new CircularTransition instance with the computed radius
        return CircularTransition(P0, P1, P2, R_star)

    @staticmethod
    def ConstrainTangentSegment(P0, P1, P2, tA, tB):
        """
        Computes the unique radius R* such that the circular transition is tangent to the constraining tangent segment.
        """
        P0 = np.array(P0)
        P1 = np.array(P1)
        P2 = np.array(P2)
        
        # Calculate PA and PB using parametric equations
        PA = P1 + tA * (P0 - P1)
        PB = P1 + tB * (P2 - P1)
        
        # Calculate the area of the triangle PA, P1, PB
        def area(PA, PB, PC):
            return 0.5 * np.abs(PA[0] * (PB[1] - PC[1]) + PB[0] * (PC[1] - PA[1]) + PC[0] * (PA[1] - PB[1]))
        
        triangle_area_PA_P1_PB = area(PA, P1, PB)
        
        # Calculate the perimeter of the triangle PA, P1, PB
        side_PA_P1 = np.linalg.norm(PA - P1)
        side_PB_P1 = np.linalg.norm(PB - P1)
        side_PA_PB = np.linalg.norm(PA - PB)
        triangle_perimeter_PA_P1_PB = side_PA_P1 + side_PB_P1 + side_PA_PB
        
        # Calculate the semiperimeter of the triangle
        triangle_semiperimeter_PA_P1_PB = triangle_perimeter_PA_P1_PB / 2
        
        # Calculate the unique radius R*
        R_star = triangle_area_PA_P1_PB / (triangle_semiperimeter_PA_P1_PB - side_PA_PB)
        print(f"R_star: {R_star}")
        
        # Return a new CircularTransition instance with the computed radius
        return CircularTransition(P0, P1, P2, R_star), PA, PB