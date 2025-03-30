# -*- coding: utf-8 -*-
"""
GeometricRoadDesign - Road Design Functions
===========================================
This module is part of the GeometricRoadDesign project and provides tools 
for analyzing and designing road alignments, including:

1. **Planimetric Curve Transitions**
   - Circular transitions between straight road sections.
   - Clothoid-based transitions with offset computation and center correction.

2. **Altimetric Curve Transitions (Vertical Alignments)**
   - Implementation planned.

3. **Complete Road Design Verification**
   - Future tools for checking geometric and regulatory compliance.

Author: Francisco
"""

import numpy as np
import plotly.graph_objects as go
from scipy.special import fresnel


class CircularTransition:
    def __init__(self, P0, P1, P2, R):
        self.P0 = np.array(P0)
        self.P1 = np.array(P1)
        self.P2 = np.array(P2)
        self.R = R

        if self.areCollinear(self.P0, self.P1, self.P2):
            raise ValueError("The three points are collinear. No circular transition is possible.")

        self.is_clockwise = self.check_curve_direction(self.P0, self.P1, self.P2)

        self.v1 = (self.P0 - self.P1) / np.linalg.norm(self.P1 - self.P0)
        self.v2 = (self.P2 - self.P1) / np.linalg.norm(self.P2 - self.P1)

        self.alpha_deg = self.ComputeIntersectionAngle()
        self.theta_deg = 180 - self.alpha_deg

        self.T = self.R * np.tan(np.radians(self.theta_deg) / 2)
        self.B = self.R / np.cos(np.radians(self.theta_deg) / 2)

        self.O = self.ComputeCircleCenter()
        self.T1, self.T2 = self.ComputeTangentPoints()

        self.arc_x, self.arc_y = self.ComputeArcCoordinates()

        self.clothoid_entry_global = None
        self.clothoid_exit_global = None
        self.O_star = None
        self.offset_line1 = None
        self.offset_line2 = None

        self.PrintParameters()

    @staticmethod
    def areCollinear(P0, P1, P2):
        area = 0.5 * np.linalg.norm(np.cross(P1 - P0, P2 - P0))
        return np.isclose(area, 0)

    @staticmethod
    def check_curve_direction(P0, P1, P2):
        cross_product = np.cross(P1 - P0, P2 - P1)
        return cross_product < 0

    def ComputeIntersectionAngle(self):
        dot_product = np.dot(self.v1, self.v2)
        alpha_rad = np.arccos(dot_product)
        return np.degrees(alpha_rad)

    def ComputeCircleCenter(self):
        bisector = (self.v1 + self.v2)
        bisector /= np.linalg.norm(bisector)
        return self.P1 + bisector * self.B

    def ComputeTangentPoints(self):
        T1 = self.P1 + self.v1 * self.T
        T2 = self.P1 + self.v2 * self.T
        return T1, T2

    def ComputeArcCoordinates(self):
        theta1 = np.arctan2(self.T1[1] - self.O[1], self.T1[0] - self.O[0])
        theta2 = np.arctan2(self.T2[1] - self.O[1], self.T2[0] - self.O[0])

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
    
    @staticmethod
    def FixPoint(P0, P1, P2, Pc):
        """
        Computes the unique radius R* such that the circular transition passes through the predefined point Pc.
        Returns a new CircularTransition instance with this radius.
        """
        P0 = np.array(P0)
        P1 = np.array(P1)
        P2 = np.array(P2)
        Pc = np.array(Pc)
    
        if CircularTransition.areCollinear(P0, P1, P2):
            raise ValueError("The three points are collinear. No circular transition is possible.")
    
        if not CircularTransition.isInsideTriangle(P0, P1, P2, Pc):
            raise ValueError("The point Pc is not inside the triangle formed by P0, P1, and P2.")
    
        v1 = (P0 - P1) / np.linalg.norm(P1 - P0)
        v2 = (P2 - P1) / np.linalg.norm(P2 - P1)
    
        alpha = np.arccos(np.dot(v1, v2))
        P1Pc = Pc - P1
        P1O = (v1 + v2) / np.linalg.norm(v1 + v2)
        omega = np.arccos(np.dot(P1Pc / np.linalg.norm(P1Pc), P1O))
    
        phi = np.arcsin(np.sin(omega) / np.sin(alpha / 2))
        if phi < np.pi / 2:
            phi = np.pi - phi
    
        gamma = np.pi - phi - omega
        R_star = np.linalg.norm(Pc - P1) * np.sin(omega) / np.sin(gamma)
    
        return CircularTransition(P0, P1, P2, R_star)
    
    @staticmethod
    def isInsideTriangle(P0, P1, P2, Pc):
        """
        Determines if the point Pc is inside the triangle formed by P0, P1, and P2.
        """
        def area(A, B, C):
            return 0.5 * np.abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
    
        A = area(P0, P1, P2)
        A1 = area(Pc, P1, P2)
        A2 = area(P0, Pc, P2)
        A3 = area(P0, P1, Pc)
    
        return np.isclose(A, A1 + A2 + A3)

    @staticmethod
    def ConstrainTangentSegment(P0, P1, P2, tA, tB):
        """
        Computes the unique radius R* such that the circular transition is tangent to the segment PA-PB.
    
        Parameters:
            P0, P1, P2: Control points defining the curve.
            tA, tB: Parameters for interpolation along segments P0P1 and P1P2.
    
        Returns:
            CircularTransition instance with computed R*, and points PA, PB defining the constraint segment.
        """
        P0 = np.array(P0)
        P1 = np.array(P1)
        P2 = np.array(P2)
    
        # Calculate points PA and PB along P0P1 and P1P2
        PA = P1 + tA * (P0 - P1)
        PB = P1 + tB * (P2 - P1)
    
        # Compute the triangle area PA-P1-PB
        def area(A, B, C):
            return 0.5 * np.abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
    
        triangle_area = area(PA, P1, PB)
    
        # Side lengths
        a = np.linalg.norm(PA - P1)
        b = np.linalg.norm(PB - P1)
        c = np.linalg.norm(PA - PB)
    
        s = (a + b + c) / 2  # Semi-perimeter
    
        # Unique radius R* (inradius of triangle)
        R_star = triangle_area / (s - c)
    
        # Build transition with R_star
        transition = CircularTransition(P0, P1, P2, R_star)
    
        return transition, PA, PB

    def PrintParameters(self):
        print(f"P0: {self.P0}, P1: {self.P1}, P2: {self.P2}")
        print(f"Radius (R): {self.R}")
        print(f"Intersection Angle α: {self.alpha_deg:.2f} degrees")
        print(f"Central Angle θ: {self.theta_deg:.2f} degrees")
        print(f"Tangent Length T: {self.T:.4f}")
        print(f"Bisector Length B: {self.B:.4f}")
        print(f"Circle Center (O): {self.O}")
        print(f"Tangent Points: T1 = {self.T1}, T2 = {self.T2}")
        
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
                                 name='Segment P0P1'))
    
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P2[0]], y=[self.P1[1], self.P2[1]],
                                 mode='lines', line=dict(color='black', width=2, dash='dash'),
                                 name='Segment P1P2'))
    
        # Reference segments from center to tangent points
        fig.add_trace(go.Scatter(x=[self.O[0], self.T1[0]], y=[self.O[1], self.T1[1]],
                                 mode='lines', line=dict(color='gray', width=1, dash='dot'),
                                 name='O–T1'))
        fig.add_trace(go.Scatter(x=[self.O[0], self.T2[0]], y=[self.O[1], self.T2[1]],
                                 mode='lines', line=dict(color='gray', width=1, dash='dot'),
                                 name='O–T2'))
    
        # Reference points
        ref_x = [self.P0[0], self.P1[0], self.P2[0], self.O[0], self.T1[0], self.T2[0]]
        ref_y = [self.P0[1], self.P1[1], self.P2[1], self.O[1], self.T1[1], self.T2[1]]
        labels = ['P0', 'P1', 'P2', 'O', 'T1', 'T2']
        fig.add_trace(go.Scatter(x=ref_x, y=ref_y, mode='markers+text',
                                 marker=dict(size=10, color='red'),
                                 text=labels, textposition="top center", name='Reference Points'))
    
        # Plot Pc if provided
        if Pc is not None:
            fig.add_trace(go.Scatter(
                x=[Pc[0]], y=[Pc[1]],
                mode='markers+text',
                marker=dict(size=10, color='green'),
                text=['Pc'], textposition="top center",
                name='Control Point'
            ))
    
        # Plot segment PA–PB if provided
        if PA is not None and PB is not None:
            fig.add_trace(go.Scatter(
                x=[PA[0], PB[0]], y=[PA[1], PB[1]],
                mode='lines+markers+text',
                line=dict(color='orange', width=2),
                marker=dict(size=10, color='orange'),
                text=['PA', 'PB'],
                textposition="top center",
                name='Segment PA–PB'
            ))
    
        # Plot offset lines if available
        if self.offset_line1 is not None:
            fig.add_trace(go.Scatter(x=[self.offset_line1[0][0], self.offset_line1[1][0]],
                                     y=[self.offset_line1[0][1], self.offset_line1[1][1]],
                                     mode='lines', line=dict(color='blue', dash='dot'),
                                     name='Offset Line 1'))
    
        if self.offset_line2 is not None:
            fig.add_trace(go.Scatter(x=[self.offset_line2[0][0], self.offset_line2[1][0]],
                                     y=[self.offset_line2[0][1], self.offset_line2[1][1]],
                                     mode='lines', line=dict(color='blue', dash='dot'),
                                     name='Offset Line 2'))
    
        # Plot O* if available
        if self.O_star is not None:
            fig.add_trace(go.Scatter(x=[self.O_star[0]], y=[self.O_star[1]],
                                     mode='markers+text',
                                     marker=dict(size=12, color='blue', symbol='x'),
                                     text=['O*'], textposition="top center",
                                     name='Corrected Center O*'))
    
        # Final layout
        fig.update_layout(title="Circular Transition Geometry",
                          xaxis_title="X", yaxis_title="Y",
                          width=900, height=700, template="plotly_white")
    
        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
    
        fig.show()


    @staticmethod
    def compute_delta_R(A, R, n_terms=6):
        L = A**2 / R
        tau = (L**2) / (2 * A**2)
        delta_sum = 0.0
        for i in range(1, n_terms + 1):
            sign = (-1)**(i + 1)
            numerator = 6 * tau**(2 * i - 2)
            denominator = (4 * i - 1) * np.math.factorial(2 * i)
            delta_sum += sign * (numerator / denominator)
        return (A**4 / (24 * R**3)) * delta_sum


    @staticmethod
    def delta_R_convergence(A, R, max_terms=6):
        """
        Analizza la convergenza della serie per il calcolo dello scostamento ΔR 
        tra clothoide e arco circolare, aumentando il numero di termini.
    
        Parametri:
            A (float): Parametro della clothoide.
            R (float): Raggio della curva circolare.
            max_terms (int): Numero massimo di termini da considerare (default: 6).
    
        Stampa:
            Una tabella con il valore di ΔR per ogni numero di termini da 1 a max_terms,
            evidenziando il contributo del 3° termine.
        """
        print(f"{'Termini':>7} | {'ΔR (m)':>12}")
        print("-" * 22)
    
        precedente = 0.0
        contributo_terzo = None
    
        for n in range(1, max_terms + 1):
            corrente = CircularTransition.compute_delta_R(A, R, n_terms=n)
            print(f"{n:>7} | {corrente:12.9f}")
            if n == 3:
                contributo_terzo = corrente - precedente
            precedente = corrente
    
        if contributo_terzo is not None:
            print(f"\n[INFO] Contributo del 3° termine: {contributo_terzo:.9e} m")

    def compute_clothoid_coordinates(self, A, is_entry):
        L = A / (self.R * np.pi)
        t_values = np.linspace(0, L / A, 100)
        S, C = fresnel(t_values)
        x_local = A * C
        y_local = A * S

        if self.is_clockwise:
            direction = (1, -1) if is_entry else (-1, -1)
        else:
            direction = (1, 1) if is_entry else (-1, 1)

        x_local *= direction[0]
        y_local *= direction[1]

        return np.vstack((x_local, y_local)).T
    
    def add_clothoids(self, A1, A2):
        print("\nComputing clothoid transition curves...")
    
        # Local clothoid coordinates from Fresnel integrals
        clothoid_entry_local = self.compute_clothoid_coordinates(A1, is_entry=True)
        clothoid_exit_local = self.compute_clothoid_coordinates(A2, is_entry=False)
    
        # ENTRY clothoid: local system aligned with P0-P1
        mid_P0P1 = (self.P0 + self.P1) / 2
        theta_P0P1 = np.arctan2(self.P1[1] - self.P0[1], self.P1[0] - self.P0[0])
        R_entry = np.array([
            [np.cos(theta_P0P1), -np.sin(theta_P0P1)],
            [np.sin(theta_P0P1),  np.cos(theta_P0P1)]
        ])
        self.clothoid_entry_global = (R_entry @ clothoid_entry_local.T).T + mid_P0P1
    
        # EXIT clothoid: local system aligned with P1-P2
        mid_P1P2 = (self.P1 + self.P2) / 2
        theta_P1P2 = np.arctan2(self.P2[1] - self.P1[1], self.P2[0] - self.P1[0])
        R_exit = np.array([
            [np.cos(theta_P1P2), -np.sin(theta_P1P2)],
            [np.sin(theta_P1P2),  np.cos(theta_P1P2)]
        ])
        self.clothoid_exit_global = (R_exit @ clothoid_exit_local.T).T + mid_P1P2
    
        # Compute ΔR offsets
        delta_R1 = self.compute_delta_R(A1, self.R)
        delta_R2 = self.compute_delta_R(A2, self.R)
    
        # Normals (pointing toward the center of the curve)
        if self.is_clockwise:
            normal1 = np.array([-self.v1[1], self.v1[0]])   # Right of v1
            normal2 = np.array([self.v2[1], -self.v2[0]])   # Right of v2
        else:
            normal1 = np.array([self.v1[1], -self.v1[0]])   # Left of v1
            normal2 = np.array([-self.v2[1], self.v2[0]])   # Left of v2
    
        self.normal1 = normal1
        self.normal2 = normal2
    
        # Parallel lines shifted from P1 along normals
        P_offset1 = self.P1 + normal1 * (self.R + delta_R1)
        P_offset2 = self.P1 + normal2 * (self.R + delta_R2)
    
        A_mat = np.column_stack((self.v1, -self.v2))
        b_vec = P_offset2 - P_offset1
    
        try:
            lambdas = np.linalg.solve(A_mat, b_vec)
            O_star = P_offset1 + lambdas[0] * self.v1
            self.O_star = O_star
            self.offset_line1 = (P_offset1 - 80 * self.v1, P_offset1 + 80 * self.v1)
            self.offset_line2 = (P_offset2 - 80 * self.v2, P_offset2 + 80 * self.v2)
            print(f"[INFO] Corrected center O* computed: {O_star}")
        except np.linalg.LinAlgError:
            print("[ERROR] Unable to compute corrected center O*. Lines are parallel.")
            self.O_star = None
            self.offset_line1 = None
            self.offset_line2 = None
            
        # === CORRECTED ARC ===
        if self.O_star is not None:
            # Clothoid lengths and characteristic angles
            L1 = A1**2 / self.R
            L2 = A2**2 / self.R
            tau1 = L1 / (2 * self.R)
            tau2 = L2 / (2 * self.R)
        
            # Use precomputed inward normal vectors
            angle1 = np.arctan2(self.normal1[1], self.normal1[0])
            angle2 = np.arctan2(self.normal2[1], self.normal2[0])
        
            # Define start and end angles
            theta_start = angle1 - tau1  # ingresso
            theta_end   = angle2 + tau2  # uscita
        
            # Ensure correct plotting direction based on curve orientation
            if self.is_clockwise:
                theta_vals = np.linspace(theta_end, theta_start, 100)  # clockwise: decreasing
            else:
                theta_vals = np.linspace(theta_start, theta_end, 100)  # counterclockwise: increasing
        
            # Save corrected arc
            self.arc_x_star = self.O_star[0] + self.R * np.cos(theta_vals)
            self.arc_y_star = self.O_star[1] + self.R * np.sin(theta_vals)
        
            theta_star = abs(theta_end - theta_start)
            print(f"[INFO] Corrected arc generated with angle θ* = {np.degrees(theta_star):.2f}°")
        else:
            self.arc_x_star = None
            self.arc_y_star = None


        print("Clothoid computation completed.\n")


    def plotClothoid(self):
        if self.clothoid_entry_global is None or self.clothoid_exit_global is None:
            raise ValueError("Clothoids not computed. Call add_clothoids(A1, A2) first.")
    
        fig = go.Figure()
    
        # Base geometry: arc and segments
        fig.add_trace(go.Scatter(x=self.arc_x, y=self.arc_y, mode='lines', name='Arc', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=[self.P0[0], self.P1[0]], y=[self.P0[1], self.P1[1]],
                                 mode='lines', name='Segment P0-P1', line=dict(color='black', dash='dash')))
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P2[0]], y=[self.P1[1], self.P2[1]],
                                 mode='lines', name='Segment P1-P2', line=dict(color='black', dash='dash')))
    
        # Clothoid curves
        fig.add_trace(go.Scatter(x=self.clothoid_entry_global[:, 0], y=self.clothoid_entry_global[:, 1],
                                 mode='lines', name='Entry Clothoid', line=dict(color='red', dash='dot')))
        fig.add_trace(go.Scatter(x=self.clothoid_exit_global[:, 0], y=self.clothoid_exit_global[:, 1],
                                 mode='lines', name='Exit Clothoid', line=dict(color='red', dash='dot')))
    
        # Reference points
        labels = ['P0', 'P1', 'P2', 'O', 'T1', 'T2']
        coords = [self.P0, self.P1, self.P2, self.O, self.T1, self.T2]
        fig.add_trace(go.Scatter(
            x=[pt[0] for pt in coords],
            y=[pt[1] for pt in coords],
            mode='markers+text',
            marker=dict(size=8, color='darkred'),
            text=labels,
            textposition='top center',
            name='Key Points'
        ))
    
        # Corrected center O*
        if self.O_star is not None:
            fig.add_trace(go.Scatter(x=[self.O_star[0]], y=[self.O_star[1]],
                                     mode='markers+text',
                                     marker=dict(size=12, color='blue', symbol='x'),
                                     text=['O*'], textposition="top center", name='Corrected Center'))
    
        # Offset lines
        if self.offset_line1 is not None:
            fig.add_trace(go.Scatter(x=[self.offset_line1[0][0], self.offset_line1[1][0]],
                                     y=[self.offset_line1[0][1], self.offset_line1[1][1]],
                                     mode='lines', line=dict(color='blue', dash='dot'), name='Offset Line 1'))
        if self.offset_line2 is not None:
            fig.add_trace(go.Scatter(x=[self.offset_line2[0][0], self.offset_line2[1][0]],
                                     y=[self.offset_line2[0][1], self.offset_line2[1][1]],
                                     mode='lines', line=dict(color='blue', dash='dot'), name='Offset Line 2'))
    
        # Tangent vectors
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.v1[0]],
                                 y=[self.P1[1], self.P1[1] + 10 * self.v1[1]],
                                 mode='lines+markers', name='v1', line=dict(color='orange')))
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.v2[0]],
                                 y=[self.P1[1], self.P1[1] + 10 * self.v2[1]],
                                 mode='lines+markers', name='v2', line=dict(color='orange')))
    
        # Normal vectors
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.normal1[0]],
                                 y=[self.P1[1], self.P1[1] + 10 * self.normal1[1]],
                                 mode='lines+markers', name='n1', line=dict(color='green')))
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.normal2[0]],
                                 y=[self.P1[1], self.P1[1] + 10 * self.normal2[1]],
                                 mode='lines+markers', name='n2', line=dict(color='green')))
        # Corrected arc (computed from O*)
        if self.arc_x_star is not None and self.arc_y_star is not None:
            fig.add_trace(go.Scatter(
                x=self.arc_x_star,
                y=self.arc_y_star,
                mode='lines',
                name='Corrected Arc (O*)',
                line=dict(color='purple', dash='solid', width=3)
            ))
    
        # Final layout
        fig.update_layout(title="Clothoid Transition with Corrected Center and Vectors",
                          xaxis_title="X", yaxis_title="Y", width=950, height=700, template="plotly_white")
        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        fig.show()
