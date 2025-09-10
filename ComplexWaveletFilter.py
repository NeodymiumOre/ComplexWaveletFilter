#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:28:14 2024

@author: joshuamarcus
"""

import os
import numpy as np
from PIL import Image
import dtcwt
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
from matplotlib import colors
import math
import tifffile as tiff
import time
import tkinter as tk
from tkinter import simpledialog, filedialog
from scipy.optimize import minimize, NonlinearConstraint
from pathlib import Path

# define parameters
N_FRAMES = 24
LAMBDA_FIRST_BIN = 405
LAMBDA_LAST_BIN = 635
LAMBDA_MIN = 400
LAMBDA_MAX = 640
IMAGE_WIDTH = 256
IMAGE_HEIGHT = 256
N_HARMONIC = 1

# COMPLEX WAVELET FILTER

# Function to select the input folder containing the G, S, and intensity .tif files
def select_input_folder():
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    input_folder = filedialog.askdirectory(title="Select the input folder containing G.tif, S.tif, and intensity.tif files")
    if input_folder:
        print(f"Selected input folder: {input_folder}")
        return input_folder
    else:
        raise ValueError("No input folder selected.")

# Function to prompt user for harmonic (H), tau (target tau value), and flevel
def get_user_inputs():
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    H = float(simpledialog.askstring("Input", "Enter Harmonic (H):"))
    tau = float(simpledialog.askstring("Input", "Enter Target Tau Value (τ):"))
    flevel = int(simpledialog.askstring("Input", "Enter Number of Filtering Levels (flevel):"))
    return H, tau, flevel

# Function to calculate Gc and Sc
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
    print(f"Calculated Gc: {Gc}, Sc: {Sc}")
    return Gc, Sc

# File handling functions
def load_and_process_image(file_path):
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    image = Image.open(file_path)  # Use PIL to load the image
    image_array = np.array(image)  # Convert to NumPy array
    return np.nan_to_num(image_array)  # Handle NaN values by replacing them with zeros


# Complex wavelet filter calculation functions
def anscombe_transform(data):
    return 2 * np.sqrt(data + (3/8))

def perform_dtcwt_transform(data, N):
    transform = dtcwt.Transform2d(biort='legall', qshift='qshift_a')
    transformed_data = transform.forward(data, nlevels=N, include_scale=False)
    return transformed_data, transform

def calculate_median_values(transformed_data):
    median_values = []
    for level in range(len(transformed_data.highpasses)):
        highpasses = transformed_data.highpasses[level if level else 0]
        for band in range(6):
            coeffs = highpasses[:, :, band]
            median_absolute = np.median(np.abs(coeffs.flatten()))
            median_values.append(median_absolute)
    return np.mean(median_values)

def calculate_local_noise_variance(transformed_data, N):
    sigma_n_squared_matrices = []

    def local_noise_variance(coeffs, N):
        sigma_n_squared = np.zeros_like(coeffs, dtype=float)
        height, width = coeffs.shape
        for x in range(width):
            for y in range(height):
                x_min, x_max = max(0, x - N), min(width, x + N + 1)
                y_min, y_max = max(0, y - N), min(height, y + N + 1)
                window = coeffs[y_min:y_max, x_min:x_max]
                local_variance = np.mean(np.abs(window)**2)
                sigma_n_squared[y, x] = local_variance
        return sigma_n_squared

    num_levels = len(transformed_data.highpasses)
    num_bands = 6

    for level in range(num_levels):
        highpasses = transformed_data.highpasses[level]
        for band in range(num_bands):
            coeffs = highpasses[:, :, band]
            sigma_n_squared = local_noise_variance(coeffs, N)
            sigma_n_squared_matrices.append((level, band, sigma_n_squared))

    return sigma_n_squared_matrices

def compute_phi_prime(mandrill_t, sigma_g_squared, sigma_n_squared_matrices):
    updated_coefficients = []

    max_level = len(mandrill_t.highpasses) - 1
    local_term = np.sqrt(3) * np.sqrt(sigma_g_squared)

    for level in range(max_level):
        highpasses_l = mandrill_t.highpasses[level]
        highpasses_l_plus_1 = mandrill_t.highpasses[level + 1]
        level_coefficients = []

        for band in range(6):
            phi_l_b = highpasses_l[:, :, band]
            phi_l_plus_1_b = highpasses_l_plus_1[:, :, band]

            _, _, sigma_n_squared = sigma_n_squared_matrices[level * 6 + band]
            phi_prime = np.zeros_like(phi_l_b, dtype=complex)

            downsample_factor = phi_l_b.shape[0] // sigma_n_squared.shape[0]

            for x in range(phi_l_b.shape[1]):
                for y in range(phi_l_b.shape[0]):
                    x_half = x // 2
                    y_half = y // 2

                    x_downsampled = x // downsample_factor
                    y_downsampled = y // downsample_factor

                    phi_squared_sum = np.abs(phi_l_b[y, x])**2 + np.abs(phi_l_plus_1_b[y_half, x_half])**2

                    if sigma_n_squared[y_downsampled, x_downsampled] > 0 and phi_squared_sum > 0:
                        denominator = np.sqrt(phi_squared_sum + local_term)
                        factor = 1 - local_term / denominator
                    else:
                        factor = 0

                    factor = max(factor, 0)
                    phi_prime[y, x] = factor * phi_l_b[y, x]

            level_coefficients.append(phi_prime)
        updated_coefficients.append(level_coefficients)

    return updated_coefficients

def update_coefficients(mandrill_t, phi_prime_matrices):
    for level, level_matrices in enumerate(phi_prime_matrices):
        for band, phi_prime in enumerate(level_matrices):
            mandrill_t.highpasses[level][:, :, band] = phi_prime

def perform_inverse_dtcwt_transform(transformed_data):
    transform = dtcwt.Transform2d(biort='legall', qshift='qshift_a')
    return transform.inverse(transformed_data)
    
def reverse_anscombe_transform(y):
    y = np.asarray(y, dtype=np.float64)
    inverse = (
        (y**2 / 4) +
        (np.sqrt(3/2) * (1/y) / 4) -
        (11 / (8 * y**2)) +
        (np.sqrt(5/2) * (1/y**3) / 8) -
        (1 / (8 * y**4))
    )
    return inverse

# GMM filtering functions
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

# Function #4: check if points are inside the ellipse
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

def check_either_value_greater_than_zero(list1, list2):
    results = [x > 0 or y > 0 for x, y in zip(list1, list2)]
    return results

def convert_list_to_array_with_dimensions(lst, rows, columns):
    array = np.array(lst)
    array_with_dimensions = array.reshape(rows, columns)
    return array_with_dimensions

# Complex wavelet filter batch processing function
def process_files(file_paths, G_combined, S_combined, I_combined, G_unfil, S_unfil, Intensity):
    #G_unfil = load_and_process_image(file_paths["G"])
    #S_unfil = load_and_process_image(file_paths["S"])
    #Intensity = load_and_process_image(file_paths["intensity"])

    if G_unfil is None or S_unfil is None or Intensity is None:
        print("One or more files could not be loaded. Skipping this replicate.")
        return G_combined, S_combined, I_combined
    
    # Compute Fourier coefficients
    Freal_rescale = G_unfil * Intensity
    Fimag_rescale = S_unfil * Intensity

    print(np.min(Freal_rescale), np.min(Fimag_rescale))
    print(np.max(Freal_rescale), np.max(Fimag_rescale))
    #exit()

    # Freal transformations and filtering
    Freal_ans = anscombe_transform(Freal_rescale)
    Freal_transformed, Freal_transformed_object = perform_dtcwt_transform(Freal_ans, flevel)
    median_values = calculate_median_values(Freal_transformed)
    sigma_g_squared = median_values / 0.6745
    sigma_n_squared = calculate_local_noise_variance(Freal_transformed, flevel)
    phi_prime = compute_phi_prime(Freal_transformed, sigma_g_squared, sigma_n_squared)
    update_coefficients(Freal_transformed, phi_prime)
    Freal_reconstructed_filtered = perform_inverse_dtcwt_transform(Freal_transformed)
    Freal_filtered = reverse_anscombe_transform(Freal_reconstructed_filtered)

    # Fimag transformations and filtering
    Fimag_ans = anscombe_transform(Fimag_rescale)
    Fimag_transformed, Fimag_transformed_object = perform_dtcwt_transform(Fimag_ans, flevel)
    median_values = calculate_median_values(Fimag_transformed)
    sigma_g_squared = median_values / 0.6745
    sigma_n_squared = calculate_local_noise_variance(Fimag_transformed, flevel)
    phi_prime = compute_phi_prime(Fimag_transformed, sigma_g_squared, sigma_n_squared)
    update_coefficients(Fimag_transformed, phi_prime)
    Fimag_reconstructed_filtered = perform_inverse_dtcwt_transform(Fimag_transformed)
    Fimag_filtered = reverse_anscombe_transform(Fimag_reconstructed_filtered)

    # Intensity transformations and filtering
    Intensity_ans = anscombe_transform(Intensity)
    Intensity_transformed, Intensity_transformed_object = perform_dtcwt_transform(Intensity_ans, flevel)
    median_values = calculate_median_values(Intensity_transformed)
    sigma_g_squared = median_values / 0.6745
    sigma_n_squared = calculate_local_noise_variance(Intensity_transformed, flevel)
    phi_prime = compute_phi_prime(Intensity_transformed, sigma_g_squared, sigma_n_squared)
    update_coefficients(Intensity_transformed, phi_prime)
    Intensity_reconstructed_filtered = perform_inverse_dtcwt_transform(Intensity_transformed)
    Intensity_filtered = reverse_anscombe_transform(Intensity_reconstructed_filtered)

    G_wavelet_filtered = Freal_filtered / Intensity_filtered
    S_wavelet_filtered = Fimag_filtered / Intensity_filtered

    threshold = 0
    
    G_array = np.nan_to_num(G_wavelet_filtered)
    S_array = np.nan_to_num(S_wavelet_filtered)
    I_array = np.nan_to_num(Intensity)

    threshold_array = I_array > threshold
    G001_array = np.clip(G_array * threshold_array, -10.0, 10.0)
    S001_array = np.clip(S_array * threshold_array, -10.0, 10.0)

    # Append to combined arrays
    if G_combined.size == 0:
        G_combined = G001_array
    else:
        G_combined = np.hstack((G_combined, G001_array))

    if S_combined.size == 0:
        S_combined = S001_array
    else:
        S_combined = np.hstack((S_combined, S001_array))

    if I_combined.size == 0:
        I_combined = I_array
    else:
        I_combined = np.hstack((I_combined, I_array))
    
    return G_combined, S_combined, I_combined

   

# Complex wavelet filter batch processing function
def process_unfil_files(file_paths, G_combined, S_combined, I_combined):
    G_unfil = tiff.imread(file_paths["G"])
    S_unfil = tiff.imread(file_paths["S"])
    Intensity = tiff.imread(file_paths["intensity"])

    if G_unfil is None or S_unfil is None or Intensity is None:
        print("One or more files could not be loaded. Skipping this replicate.")
        return G_combined, S_combined, I_combined

    threshold = 0
    
    G_array = np.nan_to_num(G_unfil)
    S_array = np.nan_to_num(S_unfil)
    I_array = np.nan_to_num(Intensity)

    threshold_array = I_array > threshold
    G001_array = np.clip(G_array * threshold_array, -0.1, 1.1)
    S001_array = np.clip(S_array * threshold_array, -0.1, 1.1)
    
    # Append to combined arrays
    if G_combined.size == 0:
        G_combined = G001_array
    else:
        G_combined = np.hstack((G_combined, G001_array))

    if S_combined.size == 0:
        S_combined = S001_array
    else:
        S_combined = np.hstack((S_combined, S001_array))

    if I_combined.size == 0:
        I_combined = I_array
    else:
        I_combined = np.hstack((I_combined, I_array))
    
    return G_combined, S_combined, I_combined
    
# Phasor and lifetime calculation from batch processed data
def plot_combined_data(G_combined, S_combined, I_combined, phasor_output_path):
    G_combined_flat = G_combined.ravel()
    S_combined_flat = S_combined.ravel()
    I_combined_flat = I_combined.ravel().astype(int)

    x_scale = [-1.005, 1.005]
    y_scale = [-1, 1]

    G001_weighted = np.repeat(G_combined_flat, I_combined_flat)
    S001_weighted = np.repeat(S_combined_flat, I_combined_flat)
    G001_weighted = np.nan_to_num(G001_weighted)
    S001_weighted = np.nan_to_num(S001_weighted)

    x = np.linspace(0, 1.0, 100)
    y = np.linspace(0, 1.0, 100)
    X, Y = np.meshgrid(x, y)
    F = (X**2 + Y**2 - X)

    iqr_x = np.percentile(G001_weighted, 75) - np.percentile(G001_weighted, 25)
    bin_width_x = 2 * iqr_x * (len(G001_weighted) ** (-1/3))
    bin_width_x = np.nan_to_num(bin_width_x)

    iqr_y = np.percentile(S001_weighted, 75) - np.percentile(S001_weighted, 25)
    bin_width_y = 2 * iqr_y * (len(S001_weighted) ** (-1/3))
    bin_width_y = np.nan_to_num(bin_width_y)

    num_bins_x_G001 = int(np.ceil((np.max(G001_weighted) - np.min(G001_weighted)) / bin_width_x)) // 2
    num_bins_y_G001 = int(np.ceil((np.max(S001_weighted) - np.min(S001_weighted)) / bin_width_y)) // 2

    hist_vals, _, _ = np.histogram2d(G_combined_flat, S_combined_flat, bins=(num_bins_x_G001, num_bins_y_G001), weights=I_combined_flat)
    vmax = hist_vals.max()
    vmin = hist_vals.min()

    fig, ax = plt.subplots(figsize=(8, 6))
    h = ax.hist2d(G_combined_flat, S_combined_flat, bins=(num_bins_x_G001, num_bins_y_G001), weights=I_combined_flat, cmap='nipy_spectral', norm=colors.SymLogNorm(linthresh=50, linscale=1, vmax=vmax, vmin=vmin), zorder=1, cmin=0.01)
    ax.set_facecolor('white')
    ax.set_xlabel('\n$G$')
    ax.set_ylabel('$S$\n')
    ax.set_xlim(x_scale)
    ax.set_ylim(y_scale)
    ax.contour(X, Y, F, [0], colors='black', linewidths=1, zorder=2)

    near_zero = 0.1
    cbar = fig.colorbar(h[3], ax=ax, format=LogFormatter(10, labelOnlyBase=True))
    ticks = [near_zero] + [10**i for i in range(1, int(np.log10(vmax)) + 1)]
    cbar.set_ticks(ticks)
    tick_labels = ['0'] + [f'$10^{i}$' for i in range(1, int(np.log10(vmax)) + 1)]
    cbar.set_ticklabels(tick_labels)
    cbar.set_label('Frequency')

    fig.tight_layout()
    fig.savefig(phasor_output_path, format='png', dpi=300)
    plt.show()
    
def calculate_and_plot_lifetime(G_combined, S_combined):
    """
    Calculate the fluorescence lifetime from combined G and S components and plot the result.

    Parameters:
    G_combined (numpy.ndarray): Combined G component array
    S_combined (numpy.ndarray): Combined S component array

    Returns:
    numpy.ndarray: Calculated fluorescence lifetime array
    """
    # Calculate phase angle theta
    SoverG = S_combined / G_combined
    theta1 = np.arctan(SoverG)

    # Calculate tangent of phase angle
    tantheta = np.tan(theta1)

    # 78mHz = period of laser - needs to be converted to nanoseconds (value below)
    ns = 12.820512820513

    # Calculate angular frequency (w) from period in nanoseconds (ns)
    w = (2 * np.pi) / ns

    # Fluorescence lifetime of image or first ROI
    T = tantheta / w

    # Plot fluorescence lifetime with colorbar (T)
    T[np.isnan(T)] = 0
    plt.imshow(T)
    plt.colorbar()
    plt.show()

    return T

def import_image_stack(source_dir):
    """Return np.array of size (IMAGE_HEIGHT, IMAGE_WIDTH, N_CHANNELS)
       containing pixels from txt files. One channel per file. 
       Order of channels is determined after sorting files by name."""

    # create list of txt files sorted by name
    txt_files = sorted([file for file in source_dir.iterdir() if file.suffix == ".txt"])
    print(*txt_files)

    # tensor for storing imported data
    data_tensor = np.zeros((256, 256, 24))

    # extract data points from each file
    for channel, file in enumerate(txt_files):
        channel_data = file.read_text().splitlines()
        # convert data from text lines to numpy.array and append it to data_tensor
        channel_data = np.array([line.strip().split(" ") for line in channel_data], dtype=int)
        data_tensor[:, :, channel] = channel_data

    return data_tensor

def convert_stack_into_phasor(stack_tensor):
    
    print("Converting stack into phasor...")

    gs_tensor = np.zeros((stack_tensor.shape[0], stack_tensor.shape[1], 3))
    for x in range(stack_tensor.shape[0]):
        for y in range(stack_tensor.shape[1]):
            gs_tensor[x, y, :] = transform_pixel_into_gs_form(stack_tensor[x, y, :])
    
    print(f"Converted stack of shape {stack_tensor.shape} into phasor of shape {gs_tensor.shape}")

    return gs_tensor


def transform_pixel_into_gs_form(data_vector, threshold=0):
    """Return np.array(g, s) for given pixel"""

    # calculate wavelength grid (create 2x more points to represent centers of each spectral range)
    l_vector = np.linspace(LAMBDA_MIN, LAMBDA_MAX, 2*N_FRAMES+1)
    # leave only centers of spectral ranges in l_vector (take every 2nd point starting from index=1)
    l_vector = l_vector[1::2]
    # calculate sum of intensity over all channels for selected pixel
    i_sum = sum(data_vector)
    # calculate total wavelength range
    L = LAMBDA_MAX - LAMBDA_MIN

    # filtering out 0s and noise values
    if i_sum <= threshold:
        return np.zeros(3)
    
    # calculate vectors for trigonometric functions
    arg = 2*np.pi*N_HARMONIC*(l_vector-l_vector[0])/L
    sinus = np.sin(arg)
    cosinus = np.cos(arg)

    # g and s should be divided by i_sum, but that will happen after filtering them first
    # calculate g component
    g = np.sum(data_vector*cosinus)/i_sum
    # calculate s component
    s = np.sum(data_vector*sinus)/i_sum

    # return g and s components as an array
    return np.array([g, s, i_sum])

def scatter_phasor(data, file_path):

    setup_figure()

    # plot data points
    plt.scatter(data[:, :, 0], data[:, :, 1], s=1, alpha=0.5)
    
    # setup plotting parameters
    setup_polar_plot()

    plt.savefig(file_path, bbox_inches='tight')

def setup_figure(figsize=(8, 5), scale=1.0):

    scaled_figsize = (figsize[0]*scale, figsize[1]*scale)
    plt.figure(figsize=scaled_figsize)
    plt.gca().set_aspect('equal', adjustable='box')

def setup_polar_plot(grid=True, black_bg=False):
    
    c = 'white' if black_bg else 'black'
    
    # draw a unit circle for reference
    circle = plt.Circle((0, 0), 1, color=c, fill=False)
    plt.gca().add_artist(circle)

    # setup plotting parameters
    plt.grid(grid)
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.xlabel("g component", fontsize=24)
    plt.ylabel("s component", fontsize=24)
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    
    if black_bg:
        # change background color
        plt.gca().set_facecolor('black')


# Main execution script
if __name__ == "__main__":
    # Record the start time
    start_time = time.time()

    # Select the input folder (which should contain G.tif, S.tif, and intensity.tif)
    #input_folder = select_input_folder()
    input_folder = None

    # Automatically set output directory to the current script location
    output_base_directory = os.path.dirname(os.path.abspath(__file__))
    print(f"Output directory is set to: {output_base_directory}")

    # Get harmonic, tau, and flevel from user
    #H, tau, flevel = get_user_inputs()
    H, tau, flevel = 1, None, 2

    # Import stack data from text files
    src_dir = Path("../../data/us_mouse_retina/730_nm_sample")
    stack_raw = import_image_stack(src_dir)

    # Calculate phsasor
    phasor = convert_stack_into_phasor(stack_raw)
    print(np.min(phasor[:,:,0]), np.min(phasor[:,:,1]))
    print(np.max(phasor[:,:,0]), np.max(phasor[:,:,1]))
    #exit()

    G_raw = phasor[:, :, 0] + 1
    S_raw = phasor[:, :, 1] + 1
    I_raw = phasor[:, :, 2]
    print(np.min(G_raw), np.min(S_raw))
    print(np.max(G_raw), np.max(S_raw))

    # Automatically calculate Gc and Sc
    #Gc, Sc = calculate_g_and_s(H, tau)

    # Define file paths for G, S, and intensity files in the input folder
    #g_tif = os.path.join(input_folder, "G.tif")
    #s_tif = os.path.join(input_folder, "S.tif")
    #intensity_tif = os.path.join(input_folder, "intensity.tif")

#    # Check if the required files exist
#    if not all(map(os.path.exists, [g_tif, s_tif, intensity_tif])):
#        raise FileNotFoundError("One or more of the required .tif files (G.tif, S.tif, intensity.tif) are missing in the input folder.")

    # Initialize arrays (no need to combine datasets, so these are just placeholders)
    G_combined = np.array([]).reshape(0, 0)
    S_combined = np.array([]).reshape(0, 0)
    I_combined = np.array([]).reshape(0, 0)
    G_combined_unfil = np.array([]).reshape(0, 0)
    S_combined_unfil = np.array([]).reshape(0, 0)
    I_combined_unfil = np.array([]).reshape(0, 0)

#    # Process the G, S, and intensity .tif files
#    file_paths = {
#        "G": g_tif,
#        "S": s_tif,
#        "intensity": intensity_tif
#    }
    file_paths = []
    G_combined, S_combined, I_combined = process_files(file_paths, G_combined, S_combined, I_combined, G_raw, S_raw, I_raw)
    #G_combined_unfil, S_combined_unfil, I_combined_unfil  = process_unfil_files(file_paths, G_combined_unfil, S_combined_unfil, I_combined_unfil)
    G_combined = G_combined - 1
    S_combined = S_combined - 1

    # Define output file names based on the provided flevel and processing conditions
    phasor_title = f'phasor_CWFlevels={flevel}.png'
    npz_file_name = f'dataset_CWFlevels={flevel}.npz'
    tiff_file_name = f'lifetime_CWFlevels={flevel}.tiff'
    tiff_file_name_unfil = f'lifetime_unfiltered.tiff'

    # Define output directories for different types of data
    phasors_dir = os.path.join(output_base_directory, "phasors")
    datasets_dir = os.path.join(output_base_directory, "datasets")
    lifetime_images_dir = os.path.join(output_base_directory, "lifetime_images")
    lifetime_images_unfil_dir = os.path.join(output_base_directory, "lifetime_images_unfiltered")

    # Ensure directories exist
    os.makedirs(phasors_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    os.makedirs(lifetime_images_dir, exist_ok=True)
    os.makedirs(lifetime_images_unfil_dir, exist_ok=True)

    # Plot combined data and save results
    phasor_path = os.path.join(phasors_dir, phasor_title)
    plot_combined_data(G_combined, S_combined, I_combined, phasor_path)

    # Calculate lifetime images and save results
    T = calculate_and_plot_lifetime(G_combined, S_combined)
    tiff_path = os.path.join(lifetime_images_dir, tiff_file_name)
    tiff.imwrite(tiff_path, T)
    #phasor_fil = np.array([G_combined, S_combined])
    phasor_fil_ext = [G_combined.tolist(), S_combined.tolist(), I_combined.tolist()]
    from copy import deepcopy
    phasor_fil = deepcopy(phasor)
    phasor_fil[:, :, 0] = G_combined
    phasor_fil[:, :, 1] = S_combined
    scatter_phasor(phasor_fil, "phasor.png")
    import json
    with open("data.json", "w", encoding='utf-8') as f:
        json.dump(phasor_fil_ext, f, indent=4)

    #T_unfil = calculate_and_plot_lifetime(G_combined_unfil, S_combined_unfil)
    tiff_path_unfil = os.path.join(lifetime_images_unfil_dir, tiff_file_name_unfil)
    tiff.imwrite(tiff_path_unfil, T_unfil)

    # Save the combined data in an .npz file
    npz_path = os.path.join(datasets_dir, npz_file_name)
    np.savez(npz_path, G=G_combined, S=S_combined, A=I_combined, T=T)

    print(f"Phasor plot saved as {phasor_path}")
    print(f"NPZ dataset saved as {npz_path}")
    print(f"Lifetime image saved as {tiff_path}")

    # Record the end time
    end_time = time.time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")






