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
www.franciscojmendez.com
"""

import numpy as np
import plotly.graph_objects as go

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
            print(f"\n[INFO] Contributo del 3° termine: {contributo_terzo:.9e} m\n")


    @staticmethod
    def compute_Xm(A, R, n_terms=10):
        """
        Computes the midpoint coordinate Xm of the clothoid using a series approximation.
    
        Parameters:
            A (float): Clothoid parameter.
            R (float): Radius of the circular arc.
            n_terms (int): Number of terms in the series expansion (default: 10).
    
        Returns:
            float: Estimated midpoint Xm in meters.
        """
        tau = A**2 / (2 * R**2)
        sqrt_2tau = np.sqrt(2 * tau)
        sum_series = 0.0
    
        for i in range(1, n_terms + 1):
            sign = (-1)**(i + 1)
            numerator = tau**(2 * i - 2)
            denominator = (4 * i - 3) * np.math.factorial(2 * i - 2)
            sum_series += sign * (numerator / denominator)
    
        Xm = 0.5 * A * sqrt_2tau * sum_series
        return Xm

    @staticmethod
    def Xm_convergence(A, R, max_terms=10):
        """
        Analyzes the convergence of the clothoid midpoint Xm by increasing the number of terms.
    
        Parameters:
            A (float): Clothoid parameter.
            R (float): Radius of the circular arc.
            max_terms (int): Maximum number of terms to include in the analysis (default: 10).
        """
        print(f"{'Terms':>6} | {'X_m (m)':>12}")
        print("-" * 22)
    
        previous = 0.0
        contribution_third = None
    
        for n in range(1, max_terms + 1):
            current = CircularTransition.compute_Xm(A, R, n_terms=n)
            print(f"{n:6} | {current:12.9f}")
            if n == 3:
                contribution_third = current - previous
            previous = current
    
        if contribution_third is not None:
            print(f"\n[INFO] Contribution of 3rd term: {contribution_third:.9e} m\n")

    @staticmethod
    def compute_clothoid(A, R, tau, n_terms=10):
        """
        Computes the clothoid point (X, Y) at a given tau using series approximation.
    
        Parameters:
            A (float): Clothoid parameter.
            R (float): Radius of the circular arc.
            tau (float): Normalized parameter τ = L² / (2A²), in [0, τ_f].
            n_terms (int): Number of terms in the series.
    
        Returns:
            (X, Y): Coordinates at the given τ.
        """
        sqrt_2tau = np.sqrt(2 * tau)
        sum_X = 0.0
        sum_Y = 0.0
        for i in range(1, n_terms + 1):
            sum_X += (-1)**(i + 1) * tau**(2 * i - 2) / ((4 * i - 3) * np.math.factorial(2 * i - 2))
            sum_Y += (-1)**(i + 1) * tau**(2 * i - 1) / ((4 * i - 1) * np.math.factorial(2 * i - 1))
    
        X = A * sqrt_2tau * sum_X
        Y = A * sqrt_2tau * sum_Y
        return X, Y

    @staticmethod
    def clothoid_convergence(A, R, tau, max_terms=10):
        """
        Prints the convergence table for X(τ) and Y(τ) values of the clothoid.
    
        Parameters:
            A (float): Clothoid parameter.
            R (float): Radius of curvature.
            tau (float): Evaluation point in [0, τ_f].
            max_terms (int): Maximum number of terms to include in the series.
        """
        print(f"{'Terms':>6} | {'X(τ) (m)':>12} | {'Y(τ) (m)':>12}")
        print("-" * 37)
    
        sqrt_2tau = np.sqrt(2 * tau)
    
        prev_X = 0.0
        prev_Y = 0.0
        contrib_X = None
        contrib_Y = None
    
        for n in range(1, max_terms + 1):
            sum_X = 0.0
            sum_Y = 0.0
            for i in range(1, n + 1):
                base = tau ** (2 * i - 2)
                factorial_term = np.math.factorial(2 * i - 2)
                sum_X += (-1) ** (i - 1) * base / ((4 * i - 3) * factorial_term)
                sum_Y += (-1) ** (i + 1) * tau * base * (4 * i - 3) / ((4 * i - 1) * factorial_term)
    
            X = A * sqrt_2tau * sum_X
            Y = A * sqrt_2tau * sum_Y
            print(f"{n:6} | {X:12.6f} | {Y:12.6f}")
    
            if n == 3:
                contrib_X = X - prev_X
                contrib_Y = Y - prev_Y
    
            prev_X = X
            prev_Y = Y
    
        if contrib_X is not None and contrib_Y is not None:
            print("\n[INFO] Contribution of 3rd term:")
            print(f"       ΔX₃ ≈ {contrib_X:.9e} m")
            print(f"       ΔY₃ ≈ {contrib_Y:.9e} m\n")
    
    def add_clothoids(self, A1, A2, n_points=20):
        print("\n[INFO] Computing entry and exit clothoids...")    
    
        # === Store input parameters ===
        self.A1 = A1
        self.A2 = A2
        self.n_points = n_points
    
        # === Compute ΔR, L and tau ===
        self.delta_R1 = self.compute_delta_R(A1, self.R)
        self.delta_R2 = self.compute_delta_R(A2, self.R)
        self.L1 = A1**2 / self.R
        self.L2 = A2**2 / self.R
        self.tau1 = self.L1 / (2 * self.R)
        self.tau2 = self.L2 / (2 * self.R)
    
        # === Compute normals ===
        if self.is_clockwise:
            self.normal1 = np.array([-self.v1[1], self.v1[0]])   # Right of v1
            self.normal2 = np.array([self.v2[1], -self.v2[0]])   # Right of v2
        else:
            self.normal1 = np.array([self.v1[1], -self.v1[0]])   # Left of v1
            self.normal2 = np.array([-self.v2[1], self.v2[0]])   # Left of v2
    
        # === Parallel offset lines ===
        P_offset1 = self.P1 + self.normal1 * (self.R + self.delta_R1)
        P_offset2 = self.P1 + self.normal2 * (self.R + self.delta_R2)
        A_mat = np.column_stack((self.v1, -self.v2))
        b_vec = P_offset2 - P_offset1
    
        try:
            lambdas = np.linalg.solve(A_mat, b_vec)
            self.O_star = P_offset1 + lambdas[0] * self.v1
            self.offset_line1 = (self.O_star - self.R * self.v1, self.O_star + self.R * self.v1)
            self.offset_line2 = (self.O_star - self.R * self.v2, self.O_star + self.R * self.v2)
            print(f"[INFO] Corrected center O* computed: {self.O_star}")
        except np.linalg.LinAlgError:
            print("[ERROR] Unable to compute corrected center O*. Lines are parallel.")
            self.O_star = None
            self.offset_line1 = None
            self.offset_line2 = None 
    
        # === Generate corrected arc ===
        if self.O_star is not None:
            theta_rad = np.radians(self.theta_deg)
            if self.tau1 + self.tau2 > theta_rad:
                print(f"[WARNING] Sum of clothoid deflections (τ₁ + τ₂ = {np.degrees(self.tau1 + self.tau2):.2f}°) exceeds arc angle θ = {self.theta_deg:.2f}°")
    
            angle1 = np.arctan2(-self.normal1[1], -self.normal1[0])
            angle2 = np.arctan2(-self.normal2[1], -self.normal2[0])
    
            if self.is_clockwise:
                theta_start = angle1 - self.tau1
                theta_end = angle2 + self.tau2
            else:
                theta_start = angle1 + self.tau1
                theta_end = angle2 - self.tau2
    
            theta_vals = np.linspace(theta_start, theta_end, 100)
            self.arc_x_star = self.O_star[0] + self.R * np.cos(theta_vals)
            self.arc_y_star = self.O_star[1] + self.R * np.sin(theta_vals)
    
            self.theta_star = theta_end - theta_start
            print(f"[INFO] Corrected arc generated with angle θ* = {np.degrees(self.theta_star):.2f}°")
        else:
            self.arc_x_star = None
            self.arc_y_star = None
    
        # === Entry clothoid ===
        self.Xm1 = self.compute_Xm(A1, self.R)
        self.d_star1 = np.linalg.norm(self.O_star - self.P1)
        self.L_total1 = np.sqrt(self.d_star1**2 - (self.R + self.delta_R1)**2)
        
        self.T1_star = self.P1 + self.v1 * (self.L_total1 + self.Xm1)
        self.C1_star = np.array([self.arc_x_star[0], self.arc_y_star[0]])
        
        tau_vals_1 = np.linspace(0, self.tau1, self.n_points)
        clothoid_entry = []
        
        for tau in tau_vals_1:
            x, y = self.compute_clothoid(A1, self.R, tau, n_terms=5)
            pt = self.T1_star - x * self.v1 + y * self.normal1
            clothoid_entry.append(pt)        
        self.clothoid_entry_global = np.array(clothoid_entry)
        # Force exact matching of start and end points
        self.clothoid_entry_global[0] = self.T1_star      # Ensure perfect start at T1*
        self.clothoid_entry_global[-1] = self.C1_star     # Ensure perfect end at C1*
        print(f"[INFO] Entry clothoid computed with {self.n_points} points (A₁ = {A1})")
        
        # === Exit clothoid ===
        self.Xm2 = self.compute_Xm(A2, self.R)
        self.d_star2 = np.linalg.norm(self.O_star - self.P1)
        self.L_total2 = np.sqrt(self.d_star2**2 - (self.R + self.delta_R2)**2)
        
        self.T2_star = self.P1 + self.v2 * (self.L_total2 + self.Xm2)
        self.C2_star = np.array([self.arc_x_star[-1], self.arc_y_star[-1]])
        
        tau_vals_2 = np.linspace(0, self.tau2, self.n_points)
        clothoid_exit = []
        
        for tau in tau_vals_2:
            x, y = self.compute_clothoid(A2, self.R, tau, n_terms=5)
            pt = self.T2_star - x * self.v2 + y * self.normal2
            clothoid_exit.append(pt)
        self.clothoid_exit_global = np.array(clothoid_exit)
        # Reverse point order to ensure forward progression from C2* to T2*
        self.clothoid_exit_global = np.array(clothoid_exit)[::-1]
        # Force exact matching of start and end points
        self.clothoid_exit_global[0] = self.C2_star     # Ensure perfect start at C2*
        self.clothoid_exit_global[-1]= self.T2_star     # Ensure perfect end at T2*

        print(f"[INFO] Exit clothoid computed with {self.n_points} points (A₂ = {A2})")


    def plotClothoid(self):
        fig = go.Figure()
        
        # Alignment points
        # labels = ['P0', 'P1', 'P2', 'O', 'T1', 'T2']
        # coords = [self.P0, self.P1, self.P2, self.O, self.T1, self.T2]
        labels = ['P0', 'P1', 'P2']
        coords = [self.P0, self.P1, self.P2]
        fig.add_trace(go.Scatter(
            x=[pt[0] for pt in coords],
            y=[pt[1] for pt in coords],
            mode='markers+text',
            marker=dict(size=8, color='darkred'),
            text=labels,
            textposition='top center',
            name='Alignment Points'
        ))
    
        # Base geometry: arc and segments
        #fig.add_trace(go.Scatter(x=self.arc_x, y=self.arc_y, mode='lines', name='Arc', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=[self.P0[0], self.P1[0]], y=[self.P0[1], self.P1[1]],
                                 mode='lines', name='Segment P0-P1', line=dict(color='black', dash='dash')))
        fig.add_trace(go.Scatter(x=[self.P1[0], self.P2[0]], y=[self.P1[1], self.P2[1]],
                                 mode='lines', name='Segment P1-P2', line=dict(color='black', dash='dash')))
    
        # Corrected center O*
        if self.O_star is not None:
            fig.add_trace(go.Scatter(x=[self.O_star[0]], y=[self.O_star[1]],
                                     mode='markers+text',
                                     marker=dict(size=12, color='blue', symbol='x'),
                                     text=['O*'], textposition="top center", name='Corrected Center'))
        # # Offset lines
        # if self.offset_line1 is not None:
        #     fig.add_trace(go.Scatter(x=[self.offset_line1[0][0], self.offset_line1[1][0]],
        #                              y=[self.offset_line1[0][1], self.offset_line1[1][1]],
        #                              mode='lines', line=dict(color='blue', dash='dot'), name='Offset Line 1'))
        # if self.offset_line2 is not None:
        #     fig.add_trace(go.Scatter(x=[self.offset_line2[0][0], self.offset_line2[1][0]],
        #                              y=[self.offset_line2[0][1], self.offset_line2[1][1]],
        #                              mode='lines', line=dict(color='blue', dash='dot'), name='Offset Line 2'))
    
        # # Tangent vectors
        # fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.v1[0]],
        #                           y=[self.P1[1], self.P1[1] + 10 * self.v1[1]],
        #                           mode='lines+markers', name='v1', line=dict(color='orange')))
        # fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.v2[0]],
        #                           y=[self.P1[1], self.P1[1] + 10 * self.v2[1]],
        #                           mode='lines+markers', name='v2', line=dict(color='orange')))
    
        # # Normal vectors
        # fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.normal1[0]],
        #                           y=[self.P1[1], self.P1[1] + 10 * self.normal1[1]],
        #                           mode='lines+markers', name='n1', line=dict(color='green')))
        # fig.add_trace(go.Scatter(x=[self.P1[0], self.P1[0] + 10 * self.normal2[0]],
        #                           y=[self.P1[1], self.P1[1] + 10 * self.normal2[1]],
        #                           mode='lines+markers', name='n2', line=dict(color='green')))
        
        # Corrected arc (computed from O*)
        if self.arc_x_star is not None and self.arc_y_star is not None:
            fig.add_trace(go.Scatter(
                x=self.arc_x_star,
                y=self.arc_y_star,
                mode='lines',
                name='Corrected Arc (O*)',
                line=dict(color='purple', dash='solid', width=3)
            ))
        
        # T1_star marker
        if hasattr(self, "T1_star"):
            fig.add_trace(go.Scatter(
                x=[self.T1_star[0]], y=[self.T1_star[1]],
                mode='markers+text',
                marker=dict(size=10, color='blue'),
                text=['T1*'], textposition="top center",
                name='T1_star'
            ))
        # C1_star marker
        if hasattr(self, "C1_star"):
            fig.add_trace(go.Scatter(
                x=[self.C1_star[0]], y=[self.C1_star[1]],
                mode='markers+text',
                marker=dict(size=10, color='purple'),
                text=['C1*'], textposition="top center",
                name='C1_star'
            ))
        # T2_star marker
        if hasattr(self, "T2_star"):
            fig.add_trace(go.Scatter(
                x=[self.T2_star[0]], y=[self.T2_star[1]],
                mode='markers+text',
                marker=dict(size=10, color='blue'),
                text=['T2*'], textposition="top center",
                name='T2_star'
            ))
        # C2_star marker
        if hasattr(self, "C2_star"):
            fig.add_trace(go.Scatter(
                x=[self.C2_star[0]], y=[self.C2_star[1]],
                mode='markers+text',
                marker=dict(size=10, color='purple'),
                text=['C2*'], textposition="top center",
                name='C2_star'
            ))
            
        
        # Entry clothoid (clothoid_entry_global)
        if self.clothoid_entry_global is not None:
            fig.add_trace(go.Scatter(
                x=self.clothoid_entry_global[:, 0],
                y=self.clothoid_entry_global[:, 1],
                mode='lines+markers',
                line=dict(color='green', dash='dot'),
                marker=dict(size=6, color='green'),
                name='Entry Clothoid'
            ))
        # Exit clothoid (clothoid_exit_global)
        if self.clothoid_exit_global is not None:
            fig.add_trace(go.Scatter(
                x=self.clothoid_exit_global[:, 0],
                y=self.clothoid_exit_global[:, 1],
                mode='lines+markers',
                line=dict(color='orange', dash='dot'),
                marker=dict(size=6, color='orange'),
                name='Exit Clothoid'
            ))

        # Final layout
        fig.update_layout(title="Clothoid Transition with Corrected Center and Vectors",
                          xaxis_title="X", yaxis_title="Y", width=950, height=700, template="plotly_white")
        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        fig.show()
