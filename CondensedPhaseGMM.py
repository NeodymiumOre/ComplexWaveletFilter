#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 09:54:22 2024

@author: joshuamarcus
"""

import os
import numpy as np
from PIL import Image
import dtcwt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import LogFormatter
from matplotlib import colors
import math
from sklearn.mixture import GaussianMixture
from matplotlib.patches import Ellipse, Circle
import tifffile as tiff
import time
import tkinter as tk
from tkinter import simpledialog, filedialog
from scipy.optimize import minimize, NonlinearConstraint

# CONDENSED PHASE GAUSSIAN MIXTURE MODEL SEGMENTATION

def select_npz_file():
    root = tk.Tk()
    root.withdraw()
    npz_path = filedialog.askopenfilename(
        title="Select a single .npz file for processing",
        filetypes=[("NumPy files", "*.npz")]
    )
    if npz_path:
        print(f"Selected file: {npz_path}")
        return npz_path
    else:
        raise FileNotFoundError("No file selected.")


def check_either_value_greater_than_zero(list1, list2):
    results = [x > 0 or y > 0 for x, y in zip(list1, list2)]
    return results

def load_npz_file(npz_path):

    npz_file = np.load(npz_path)
    return npz_file['G'], npz_file['S'], npz_file['A'], npz_file['T']

def get_user_inputs():
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    H = float(simpledialog.askstring("Input", "Enter Harmonic (H):"))
    tau = float(simpledialog.askstring("Input", "Enter Target Tau Value (τ):"))
    cov_f = float(simpledialog.askstring("Input", "Enter multiplication factor for the ROI size:"))
    shift = float(simpledialog.askstring("Input", "Enter shift factor:"))
    radius_ref = float(simpledialog.askstring("Input", "Enter the radius for circular condensed phase filtering ROI:"))
    
    return H, tau, cov_f, shift, radius_ref


def calculate_g_and_s(H, tau):
    # Define error function for tau
    def error_function(GS, tau, H):
        G, S = GS
        M = np.sqrt(G**2 + S**2)
        theta = np.arcsin(S / M)
        ns = 12.820512820513  # Based on your example
        w = (H * np.pi) / ns
        calculated_tau = np.tan(theta) / (w * 2)
        return (calculated_tau - tau)**2

    # Constraint function to ensure solution lies on the universal circle
    def circle_constraint(GS):
        G, S = GS
        return (G - 0.5)**2 + S**2 - 0.25

    # Generate a grid of initial guesses
    initial_guesses = [(0.5 + 0.5 * np.cos(theta), 0.5 * np.sin(theta)) for theta in np.linspace(0, 2 * np.pi, 100)]

    # Constraint for the optimizer
    circle_constraint = NonlinearConstraint(circle_constraint, 0, 0)

    # Minimize error function to find optimal G and S coordinates
    result = minimize(error_function, initial_guesses[0], args=(tau, H), constraints=[circle_constraint])
    Gc, Sc = result.x
    return Gc, Sc

def threshold_and_clip_arrays(G_array, S_array, I_array, threshold=0, clip_limits=(-0.1, 1.1)):
    G_array = G_array * (I_array > threshold)
    S_array = S_array * (I_array > threshold)
    G_array = np.nan_to_num(G_array)
    S_array = np.nan_to_num(S_array)
    G_array = np.clip(G_array, *clip_limits)
    S_array = np.clip(S_array, *clip_limits)
    return G_array, S_array

def flatten_arrays(G001_array, S001_array, I_array):
    G001 = G001_array.ravel()
    S001 = S001_array.ravel()
    I = I_array.ravel().astype(int)
    return G001, S001, I

def x_y_dimensions(I_array):
    x_dim = I_array.shape[0]
    y_dim = I_array.shape[1]
    return x_dim, y_dim

def weighted_G_and_S(G001, S001, I):
    G001_weighted = np.repeat(G001, I)
    S001_weighted = np.repeat(S001, I)
    G001_weighted = np.nan_to_num(G001_weighted)
    S001_weighted = np.nan_to_num(S001_weighted)
    return G001_weighted, S001_weighted
  
def is_point_inside_circle(point, center, radius):
    distance = math.sqrt((point[0] - center[0])**2 + (point[1] - center[1])**2)
    if distance <= radius:
        return True
    else:
        return False
    
def are_points_inside_circle(points, center, radius):
    results = []
    for point in points:
        results.append(is_point_inside_circle(point, center, radius))
    return results

def is_points_inside_rotated_ellipse(center_x, center_y, semi_major_axis, semi_minor_axis, angle_degrees, points):
    # Calculate the distance between each point and the center of the ellipse
    distances = [(point[0] - center_x)**2 + (point[1] - center_y)**2 for point in points]

    # Check if the ellipse is a circle (semi-major and semi-minor axes are equal)
    is_circle = math.isclose(semi_major_axis, semi_minor_axis)

    results = []

    if is_circle:
        # If it's a circle, check if each point is inside the circle
        for distance in distances:
            results.append(distance <= semi_major_axis**2)
    else:
        # Calculate the rotation angle of the ellipse
        angle_radians = math.radians(angle_degrees)
        cos_a = math.cos(angle_radians)
        sin_a = math.sin(angle_radians)

        for i, point in enumerate(points):
            point_x, point_y = point

            # Translate the point to the ellipse's coordinate system
            translated_x = point_x - center_x
            translated_y = point_y - center_y

            # Apply the rotation transformation
            rotated_x = cos_a * translated_x + sin_a * translated_y
            rotated_y = -sin_a * translated_x + cos_a * translated_y

            # Calculate the normalized coordinates
            normalized_x = rotated_x / semi_major_axis
            normalized_y = rotated_y / semi_minor_axis

            # Check if the transformed point is inside the unrotated ellipse
            results.append(normalized_x ** 2 + normalized_y ** 2 <= 1)

    return results

def is_points_inside_circle(center_x, center_y, radius, points):
    results = []
    for point in points:
        point_x, point_y = point
        # Calculate the distance between the point and the center of the circle
        distance = (point_x - center_x)**2 + (point_y - center_y)**2
        # Check if the point is inside the circle
        results.append(distance <= radius**2)
    
    return results

def prepare_for_gmm():
    matrix_for_cond_GMM = np.reshape(results_cond, (x_dim, y_dim))
    G_for_cond_GMM = matrix_for_cond_GMM * G001_array
    S_for_cond_GMM = matrix_for_cond_GMM * S001_array
    G_GMMcond = G_for_cond_GMM.ravel()
    S_GMMcond = S_for_cond_GMM.ravel()
    G_GMMcond_weighted = np.repeat(G_GMMcond, I)
    S_GMMcond_weighted = np.repeat(S_GMMcond, I)
    data_GMMcond = np.column_stack((G_GMMcond_weighted[G_GMMcond_weighted != 0], S_GMMcond_weighted[S_GMMcond_weighted != 0]))
    return data_GMMcond

def GMM_processing():
    num_clusters = 2
    gmm_cond = GaussianMixture(n_components=num_clusters)
    gmm_cond.fit(data_GMMcond)

    # Get cluster centers, covariances, and weights
    cluster_centers_cond = gmm_cond.means_
    cov_matrices_cond = gmm_cond.covariances_

    # Plot cluster's ellipse boundaries for cond phase
    cond_cluster_center_1 = cluster_centers_cond[0]
    cond_cluster_center_2 = cluster_centers_cond[1]
    cov_matrix_1_cond = cov_matrices_cond[0]
    cov_matrix_2_cond = cov_matrices_cond[1]
    eigenvalues_1_cond, eigenvectors_1_cond = np.linalg.eigh(cov_matrix_1_cond)
    eigenvalues_2_cond, eigenvectors_2_cond = np.linalg.eigh(cov_matrix_2_cond)

    # Reassign cluster centers and corresponding eigenvalues/eigenvectors for condensed and dilute
    real_cond_cluster_center, real_cond_eigenvalues, real_cond_eigenvectors, real_dilu_cluster_center, real_dilu_eigenvalues, real_dilu_eigenvectors = reassign_good_center(
        cond_cluster_center_1, 
        cond_cluster_center_2, 
        eigenvalues_1_cond, 
        eigenvalues_2_cond, 
        eigenvectors_1_cond, 
        eigenvectors_2_cond
        )
    return real_cond_cluster_center, real_cond_eigenvalues, real_cond_eigenvectors, real_dilu_cluster_center, real_dilu_eigenvalues, real_dilu_eigenvectors

def segmentation_ROI_parameters():
    # Compute the correct angle for the ROIs from the correct eigenvector for condition
    if real_cond_eigenvectors[0, 1] > 0 and real_cond_eigenvectors[1, 1] > 0:
        angle_cond = np.arctan2(-real_cond_eigenvectors[1, 1], -real_cond_eigenvectors[0, 1])
    else:
        angle_cond = np.arctan2(real_cond_eigenvectors[1, 1], real_cond_eigenvectors[0, 1])
    angle_degrees_cond = np.degrees(angle_cond)

    # Calculate ellipse width and height for condition
    width_cond = cov_f * np.sqrt(real_cond_eigenvalues[1])
    height_cond = cov_f * np.sqrt(real_cond_eigenvalues[0])

    # Calculate the shifts for condition
    dx_cond = shift * width_cond * np.cos(angle_cond)
    dy_cond = shift * width_cond * np.sin(angle_cond)

    # Calculate the new center coordinates for condition
    center_cond_x = real_cond_cluster_center[0] + dx_cond
    center_cond_y = real_cond_cluster_center[1] + dy_cond

    # Create a new center point for condition
    ROI_center_cond = np.array([center_cond_x, center_cond_y])
    
    # Calculating G and S coordinates inside cluster ellipses for condition
    results_cluster_cond = is_points_inside_rotated_ellipse(ROI_center_cond[0], ROI_center_cond[1], width_cond, height_cond, angle_degrees_cond, points001)
    matrix_cluster_cond = np.reshape(results_cluster_cond, (x_dim, y_dim))
    
    return ROI_center_cond, width_cond, height_cond, angle_degrees_cond, results_cluster_cond, matrix_cluster_cond

def flatten_GMM():
    G_cluster_cond = matrix_cluster_cond * G001_array
    S_cluster_cond = matrix_cluster_cond * S001_array
    G_cond = G_cluster_cond.ravel()
    S_cond = S_cluster_cond.ravel()
    return G_cond, S_cond

def reassign_good_center(cluster_center_1, cluster_center_2, eigenvalues_1, eigenvalues_2, eigenvectors_1, eigenvectors_2):
    # Extract x values
    cluster_center_1_x = cluster_center_1[0]
    cluster_center_2_x = cluster_center_2[0]
    
    # Compare x values and assign the good center and its corresponding eigenvalues and eigenvectors
    if cluster_center_1_x < cluster_center_2_x:
        good_center = cluster_center_1
        good_eigenvalues = eigenvalues_1
        good_eigenvectors = eigenvectors_1
        bad_center = cluster_center_2
        bad_eigenvalues = eigenvalues_2
        bad_eigenvectors = eigenvectors_2
    else:
        good_center = cluster_center_2
        good_eigenvalues = eigenvalues_2
        good_eigenvectors = eigenvectors_2
        bad_center = cluster_center_1
        bad_eigenvalues = eigenvalues_1
        bad_eigenvectors = eigenvectors_1
    
    return good_center, good_eigenvalues, good_eigenvectors, bad_center, bad_eigenvalues, bad_eigenvectors

def plot_2d_histogram(G001, S001, I, G001_weighted, S001_weighted, ROI_center_cond, width_cond, height_cond, angle_degrees_cond):
    # 1. Generate meshgrid for contour plot
    x = np.linspace(0, 1.0, 100)
    y = np.linspace(0, 1.0, 100)
    X, Y = np.meshgrid(x, y)
    F = (X**2 + Y**2 - X)
    
    # Define axis limits
    x_scale = [-0.005, 1.005]
    y_scale = [0, 0.9]
    
    # 2. Calculate bin widths for G001 and S001
    def calculate_bin_width(data):
        iqr = np.percentile(data, 75) - np.percentile(data, 25)
        if iqr == 0:
            print("Warning: Zero interquartile range in data. Using default bin width.")
            return 0.1  # Default value
        return 2 * iqr * (len(data) ** (-1 / 3))
    
    bin_width_x = calculate_bin_width(G001_weighted)
    bin_width_y = calculate_bin_width(S001_weighted)
    
    # Determine number of bins
    num_bins_x = int(np.ceil((np.max(G001_weighted) + np.min(G001_weighted)) / bin_width_x)) // 2
    num_bins_y = int(np.ceil((np.max(S001_weighted) + np.min(S001_weighted)) / bin_width_y)) // 2
    
    # 3. Create histogram values and set colormap limits
    hist_vals, _, _ = np.histogram2d(G001, S001, bins=(num_bins_x, num_bins_y), weights=I)
    vmax = hist_vals.max()
    vmin = hist_vals.min()
    gray_cmap = ListedColormap(['#c8c8c8'])
    
    # 4. Plot histogram with colormap
    fig, ax = plt.subplots(figsize=(8, 6))
    h = ax.hist2d(G001, S001, bins=(num_bins_x, num_bins_y), weights=I, cmap='nipy_spectral',
                  norm=colors.SymLogNorm(linthresh=100, linscale=1, vmax=vmax, vmin=vmin), zorder=1, cmin=0.01)
    ax.set_facecolor('white')
    ax.set_xlabel('\n$G$')
    ax.set_ylabel('$S$\n')
    ax.set_xlim(x_scale)
    ax.set_ylim(y_scale)
    ax.contour(X, Y, F, [0], colors='black', linewidths=1, zorder=2)
    
    # 5. Add ellipse and colorbar
    ell_cond = Ellipse(ROI_center_cond, 2 * width_cond, 2 * height_cond, angle=angle_degrees_cond,
                       color='blue', fill=False, linewidth=1)
    ax.add_patch(ell_cond)
    
    # Colorbar customization
    near_zero = 0.1
    cbar = fig.colorbar(h[3], ax=ax, format=LogFormatter(10, labelOnlyBase=True))
    ticks = [near_zero] + [10**i for i in range(1, int(np.log10(vmax)) + 1)]
    cbar.set_ticks(ticks)
    tick_labels = ['0'] + [f'$10^{i}$' for i in range(1, int(np.log10(vmax)) + 1)]
    cbar.set_ticklabels(tick_labels)
    cbar.set_label('Frequency')
    
    plt.show()
    
def convert_list_to_array_with_dimensions(lst, rows, columns):
    array = np.array(lst)
    array_with_dimensions = array.reshape(rows, columns)
    return array_with_dimensions

def GMM_segmentation_mask(x_dim, y_dim, G_cond, S_cond):
    num_rows = x_dim
    num_columns = y_dim

    # CLUSTER1 processing
    new_results_cond = check_either_value_greater_than_zero(G_cond, S_cond)
    cluster_1_list = convert_list_to_array_with_dimensions(new_results_cond, num_rows, num_columns)

    # Scale the values in matrix_dilute_45min to the range of 0-255
    rgb_cluster_cond = (cluster_1_list * 255).astype(np.uint8)

    # Convert to PIL Image
    CL1 = Image.fromarray(rgb_cluster_cond)

    # Define path and save the mask image
    mask1_path = os.path.join(os.getcwd(), 'ellipse_mask.tiff')
    CL1.save(mask1_path)

if __name__ == "__main__":
    
    # Record the start time
    start_time = time.time()

    # Load selected .npz file
    npz_path = select_npz_file()
    G_array, S_array, I_array, T = load_npz_file(npz_path)

    # Get harmonic, tau, and flevel from user
    H, tau, cov_f, shift, radius_ref = get_user_inputs()
    
    # Automatically calculate Gc and Sc
    Gc, Sc = calculate_g_and_s(H, tau)
    
    # Continue with the rest of your processing as usual
    G001_array, S001_array = threshold_and_clip_arrays(G_array, S_array, I_array)
    center_cond = (Gc, Sc)

    # flatten arrays for processing and plotting
    G001, S001, I = flatten_arrays(G001_array, S001_array, I_array)
    points001 = list(zip(G001, S001))
    
    # x and y dimensions of the datasets
    x_dim, y_dim = x_y_dimensions(I_array)

    # Weight G and S coordinates by intensity for barycenter analysis
    G001_weighted, S001_weighted = weighted_G_and_S(G001, S001, I)

    # Calculate G and S coordinates inside and outside the circle for cond phase
    results_cond = are_points_inside_circle(points001, center_cond, radius_ref)
    
    data_GMMcond = prepare_for_gmm()
    
    real_cond_cluster_center, real_cond_eigenvalues, real_cond_eigenvectors, real_dilu_cluster_center, real_dilu_eigenvalues, real_dilu_eigenvectors = GMM_processing()
    
    ROI_center_cond, width_cond, height_cond, angle_degrees_cond, results_cluster_cond, matrix_cluster_cond = segmentation_ROI_parameters()
    
    G_cond, S_cond = flatten_GMM()
     
    plot_2d_histogram(G001, S001, I, G001_weighted, S001_weighted, ROI_center_cond, width_cond, height_cond, angle_degrees_cond)
    
    GMM_segmentation_mask(x_dim, y_dim, G_cond, S_cond)

    # Record the end time
    end_time = time.time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")
