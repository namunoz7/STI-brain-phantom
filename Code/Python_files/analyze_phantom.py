#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 18:49 2023
Analyze the influence of the maximum tilting angle in the STI reconstruction. The angles vary from
[50.0, 40.0, 30.0, 20.0]. We added gaussian noise with sigma=1e-3, and we used Least Squares algorithm as STI
reconstructor

@author: Nestor Munoz
"""

import torch
import nibabel as nib
import numpy as np
import sys
import os
from scipy.io import savemat

# IMG_ROOT = '../../Imagenes/Phantom_real_data/'
# IMG_ROOT = '../Phantom_real_data/'
IMG_ROOT = '../../../researchers/cristian-tejos/datasets/Phantom_real_data/'

MASK_FILE = os.path.join(IMG_ROOT, 'Masks', 'brain_mask.nii.gz')
ATLAS_NAME = os.path.join(IMG_ROOT, 'Masks', 'atlas_mask.nii.gz')
CHI_ROOT = os.path.join(IMG_ROOT, 'Phantom_tensor_2/')
RECONSTRUCTION_FOLDERS = ['COSMOS_STI', 'STI_suite', 'STI_pytorch']
CHI_NAMES = 'chi_sti_filt.nii.gz'
MMS_GT_NAME = 'mms.nii.gz'
MSA_GT_NAME = 'msa.nii.gz'

SIMULATED_FOLDER = IMG_ROOT + 'Experiment_angles/'
SIM_THETA_NAME = 'Theta_'
SIM_PSI_NAME = 'Psi_'
PHI_BASE_NAME = 'phi_'
CHI_BASE_NAME = 'chi_'
V1_BASE_NAME = 'eig_v1_'
MMS_BASE_NAME = 'mms_'
MSA_BASE_NAME = 'msa_'

ANGLES_BASE_NAME = 'angles_'
ORIENTATIONS_BASE_NAME = '_orientations.nii.gz'

MAX_PSI = [50.0, 40.0, 30.0, 20.0]
MAX_THETA = [50.0, 40.0, 30.0, 20.0]
N_ORIENTATIONS = [12]
GAMMA = 42.58
B0 = 3  # T
STD_NOISE = 8e-4

phase_scale = 2 * np.pi * GAMMA * B0
eps = sys.float_info.epsilon


def read_img(filename):
    """
    Read nifty image
    :param filename:
    :return:
    """
    nii_img = nib.load(filename)
    img = nii_img.get_fdata()
    img = torch.from_numpy(img)
    img_shape = np.array(img.shape[0:3])
    header = nii_img.header
    voxel_size = np.array(header.get_zooms()[0:3])
    fov = img_shape * voxel_size

    return img, img_shape, fov, nii_img


def fft_phase(inputs):
    """
    Gets the fourier transform of the input geometric figure images.
    :param inputs: Input image
    :return:
    """
    fft_input = torch.fft.fftn(input=inputs, dim=(0, 1, 2))
    return fft_input


def inv_fft_phase(fft_input):
    """
    Gets the inverse Fourier Transform
    :param fft_input:
    :return:
    """
    inputs = torch.fft.ifftn(fft_input, dim=(0, 1, 2))
    return inputs


def shift_fft(x, dims, mat_size):
    """
    Shift zero-frequency component to center of spectrum
    :param mat_size:
    :param x: Input image
    :param dims: Dimensions to roll
    :return:
    """
    x = x.roll((mat_size[0] // 2,
                mat_size[1] // 2,
                mat_size[2] // 2), dims)
    return x


def angles_cylinders(num_rotations, theta, psi):
    """
    Rotation angles that are used to rotate the object in the scanner
    :param num_rotations: Number of rotations that will be scanned
    :param theta: Angle of rotation in the y-axis (LR axis)
    :param psi: Angle of rotation in the x-axis (AP axis)
    """
    min_rand = -5.0
    max_rand = 5.0

    tmp = (max_rand - min_rand) * torch.rand((1, num_rotations), dtype=torch.float64) + min_rand
    # Tilt angle of main field
    if num_rotations == 6:
        vec_theta = torch.tensor([0.0, 0.0, 0.0, theta/2, theta/2, theta], dtype=torch.float64)
    else:
        vec_theta = torch.tensor([0.0, 0.0, 0.0, theta/2, theta/2, theta/2, theta/2, theta, theta, theta, theta, theta],
                                 dtype=torch.float64)
    vec_theta = vec_theta.reshape(1, num_rotations)
    vec_theta = vec_theta + tmp
    vec_theta = torch.deg2rad(vec_theta)

    tmp = (max_rand - min_rand) * torch.rand((1, num_rotations), dtype=torch.float64) + min_rand
    # Rotation angle with the z axis
    if num_rotations == 6:
        vec_psi = torch.tensor([0.0, -psi, psi, -psi/2, psi/2, 0.0], dtype=torch.float64)
    else:
        vec_psi = torch.tensor([0.0, -psi, psi, -psi/2, psi/2, -psi, psi, -psi/2, psi/2, -psi, psi, 0.0],
                               dtype=torch.float64)
    vec_psi = vec_psi.reshape(1, num_rotations)
    vec_psi += tmp
    vec_psi = torch.deg2rad(vec_psi)

    return vec_theta, vec_psi


def get_direction_field(vec_theta, vec_psi, n_orientations):
    """
    Gets the direction field vector of the multiple orientations, made by the cylinders. All the angles are in radians
    :param vec_theta: Rotation angle in the x-z plane
    :param vec_psi: Rotation angle in the y-z plane
    :param n_orientations:
    :return:
    """
    direction_field = torch.zeros(n_orientations, 3, dtype=torch.float64)
    direction_field[:, 0] = torch.sin(vec_theta)  # Hx
    direction_field[:, 1] = -torch.sin(vec_psi) * torch.cos(vec_theta)  # Hy
    direction_field[:, 2] = torch.cos(vec_psi) * torch.cos(vec_theta)  # Hz
    direction_field = direction_field.unsqueeze(-1)

    return direction_field


def gen_k_space(fov, n_orientations, mat_size):
    """
    Defines the K space
    :param mat_size:
    :param fov:
    :param n_orientations:
    :return:
    """
    kx = torch.arange(1, mat_size[0] + 1, dtype=torch.float64)
    ky = torch.arange(1, mat_size[1] + 1, dtype=torch.float64)
    kz = torch.arange(1, mat_size[2] + 1, dtype=torch.float64)

    center_x = mat_size[0] // 2 + 1
    center_y = mat_size[1] // 2 + 1
    center_z = mat_size[2] // 2 + 1
    kx = kx - center_x
    ky = ky - center_y
    kz = kz - center_z

    delta_kx = 1 / fov[0]
    delta_ky = 1 / fov[1]
    delta_kz = 1 / fov[2]

    #  Generation of k space

    kx = kx * delta_kx
    ky = ky * delta_ky
    kz = kz * delta_kz

    kxx, kyy, kzz = torch.meshgrid(kx, ky, kz)

    kxx = kxx.unsqueeze(3).repeat(1, 1, 1, n_orientations)
    kyy = kyy.unsqueeze(3).repeat(1, 1, 1, n_orientations)
    kzz = kzz.unsqueeze(3).repeat(1, 1, 1, n_orientations)

    k = torch.zeros(mat_size[0], mat_size[1], mat_size[2], n_orientations, 3, dtype=torch.float64)
    k[:, :, :, :, 0] = kxx
    k[:, :, :, :, 1] = kyy
    k[:, :, :, :, 2] = kzz
    k = k.unsqueeze(-1)

    return k, kxx, kyy, kzz


def projection_variables(vec_theta, vec_psi, n_orientations, fov, mat_size):
    """
    Generates the projection variables of the STI model (a_ii and a_ij) and construct the projection matrix
    These are only made with 12 different orientation
    :param mat_size:
    :param fov:
    :param vec_theta: vector containing the angle of deviation with respect to the main field axis
    :param vec_psi: vector containing the angle of rotation of the x-y plane.
    :param n_orientations:
    :return: A: matrix containing each one of the projection angles.
    """
    direction_field = get_direction_field(vec_theta, vec_psi, n_orientations)
    print('...... Generating k space')
    k, kxx, kyy, kzz = gen_k_space(fov, n_orientations, mat_size)
    print('...... Done')
    k2 = (kxx * kxx) + (kyy * kyy) + (kzz * kzz)
    k2[k2 == 0] = eps
    kt_h = torch.matmul(k.transpose(-2, -1), direction_field).squeeze()
    direction_field = direction_field.squeeze()
    print('...... Calculating auxiliary variables')
    # Aux variables are defined
    a_11 = ((direction_field[:, 0] * direction_field[:, 0]) / 3) - ((kt_h / k2) * kxx * direction_field[:, 0])
    a_22 = ((direction_field[:, 1] * direction_field[:, 1]) / 3) - ((kt_h / k2) * kyy * direction_field[:, 1])
    a_33 = ((direction_field[:, 2] * direction_field[:, 2]) / 3) - ((kt_h / k2) * kzz * direction_field[:, 2])
    a_12 = ((2 / 3) * direction_field[:, 0] * direction_field[:, 1]) - \
           ((kt_h / k2) * (kxx * direction_field[:, 1] + kyy * direction_field[:, 0]))
    a_13 = ((2 / 3) * direction_field[:, 0] * direction_field[:, 2]) - \
           ((kt_h / k2) * (kxx * direction_field[:, 2] + kzz * direction_field[:, 0]))
    a_23 = ((2 / 3) * direction_field[:, 2] * direction_field[:, 1]) - \
           ((kt_h / k2) * (kzz * direction_field[:, 1] + kyy * direction_field[:, 2]))
    print('...... Done')
    matrix_projection = torch.zeros(mat_size[0], mat_size[1], mat_size[2], n_orientations, 6, dtype=torch.float64)
    matrix_projection[:, :, :, :, 0] = a_11
    matrix_projection[:, :, :, :, 1] = a_12
    matrix_projection[:, :, :, :, 2] = a_13
    matrix_projection[:, :, :, :, 3] = a_22
    matrix_projection[:, :, :, :, 4] = a_23
    matrix_projection[:, :, :, :, 5] = a_33

    matrix_projection = shift_fft(matrix_projection, (0, 1, 2), mat_size)

    return matrix_projection


def get_total_phase(chi, matrix_projection):
    """
    Calculates the total phase of the tensor using the linear form of the Tensor model
    :param chi: Vectorized form of the tensor
    :param matrix_projection: Projection matrix calculated by the corresponding rotations
    :return: phase: the local phase of the image
    """
    fft_chi = fft_phase(chi)
    tmp_real = torch.matmul(matrix_projection, torch.real(fft_chi).unsqueeze(-1))
    tmp_img = torch.matmul(matrix_projection, torch.imag(fft_chi).unsqueeze(-1))
    tmp_phi = torch.cat((tmp_real, tmp_img), dim=-1)
    phi = inv_fft_phase(torch.view_as_complex(tmp_phi))
    return torch.real(phi)


def min_squares(bulk_phase, mat_projection):
    """
    Calculates the tensor by getting the phase of geometric figures.
    :param mat_projection:
    :param bulk_phase: Total phase of the image phase
    :return: The vectorized form of the tensor
    """
    mat_transpose = mat_projection.transpose(3, 4)
    b_mat = torch.matmul(mat_transpose, mat_projection)
    b_inv = b_mat.inverse()

    fft_phi = fft_phase(bulk_phase)
    tmp_real = torch.matmul(mat_transpose, torch.real(fft_phi).unsqueeze(-1))
    tmp_img = torch.matmul(mat_transpose, torch.imag(fft_phi).unsqueeze(-1))
    real_ft_chi = torch.matmul(b_inv, tmp_real).unsqueeze(-1)
    img_ft_chi = torch.matmul(b_inv, tmp_img).unsqueeze(-1)
    tmp_chi = torch.cat((real_ft_chi, img_ft_chi), dim=-1)
    chi = inv_fft_phase(torch.view_as_complex(tmp_chi))
    chi = chi.squeeze()
    return torch.real(chi)


def save_img(nii, img, name):
    """
    Save the image as a nifti file
    :param nii: Nifti variable to save
    :param img: Tensor to save
    :param name: Name of the tensor
    """
    hdr = nii.header
    new_dim = img.shape
    hdr.set_data_shape(new_dim)
    nii = nib.Nifti1Image(img, affine=None, header=hdr)
    nib.save(nii, name)
    print('...... ' + name + ' saved')


def mask_image(img, mask):
    """
    Mask the input image with the mask, to obtain only the tissue inside the brain.
    :param img:
    :param mask:
    :return:
    """
    n_orientations = img.shape[-1]
    img = torch.mul(img, mask.unsqueeze(-1).repeat([1, 1, 1, n_orientations]))
    return img


def get_eigen_decomposition(tensor):
    """
    Generates the eigen decomposition of the STI tensor
    :param tensor: Susceptibility tensor in vector form
    :return: eigen_values and eigen_vectors
    """
    mat_size = tensor.size()
    tensor = torch.stack(tensors=[tensor[:, :, :, 0], tensor[:, :, :, 1], tensor[:, :, :, 2],
                                  tensor[:, :, :, 1], tensor[:, :, :, 3], tensor[:, :, :, 4],
                                  tensor[:, :, :, 2], tensor[:, :, :, 4], tensor[:, :, :, 5]], dim=-1)
    tensor = torch.reshape(tensor, shape=(mat_size[0], mat_size[1], mat_size[2], 3, 3))
    eig_val, eig_vec = torch.linalg.eigh(tensor)
    eig_val = torch.stack([eig_val[:, :, :, 2], eig_val[:, :, :, 1], eig_val[:, :, :, 0]], dim=-1)
    v1 = eig_vec[:, :, :, :, 0]
    return eig_val, v1


def gen_mms_msa(tensor):
    """
    Calculates the mean magnetic susceptibility (MMS) and magnetic susceptibility anisotropy (MSA)
    :param tensor: Susceptibility tensor input in vectorized form
    :return: mms, msa
    """
    eigenvalues, v1 = get_eigen_decomposition(tensor)
    mms = torch.mean(eigenvalues, dim=-1)
    msa = eigenvalues[:, :, :, 0] - torch.mean(eigenvalues[:, :, :, 1:3], dim=-1)
    return v1, mms, msa


def get_angle_images(pev_1, pev_2):
    """
    Calculates the cosine of both principal eigenvectors.
    :param pev_1:
    :param pev_2:
    :return:
    """
    pev_1 = torch.unsqueeze(pev_1, dim=-1)
    pev_2 = torch.unsqueeze(pev_2, dim=-2)
    dot_prod = torch.matmul(pev_2, pev_1)
    dot_prod = torch.abs(torch.squeeze(dot_prod))
    return torch.squeeze(dot_prod)


def check_root(root_file):
    """
    Check if the root file exists. If it does not, it creates the root file
    :param root_file:
    :return:
    """
    if not os.path.exists(root_file):
        os.mkdir(root_file)
        print("Directory " + root_file + " Created ")
    else:
        print("Directory " + root_file + " already exists")


def main():
    chi_gt_file = os.path.join(CHI_ROOT, CHI_NAMES)
    chi_gt, _, _, nii = read_img(chi_gt_file)
    v1_gt, mms_gt, msa_gt = gen_mms_msa(chi_gt)

    white_matter, _, _, _ = read_img(ATLAS_NAME)
    white_matter = white_matter == 3
    for n_theta in MAX_THETA:
        for n_psi in MAX_PSI:
            for orientation in N_ORIENTATIONS:
                print('=======================================')
                actual_theta_name = SIM_THETA_NAME + str(n_theta)
                actual_psi_name = SIM_PSI_NAME + str(n_psi)
                actual_theta_folder = os.path.join(SIMULATED_FOLDER, actual_theta_name)
                actual_angles_root = os.path.join(SIMULATED_FOLDER, actual_theta_name, actual_psi_name)
                phi_name = PHI_BASE_NAME + str(orientation) + ORIENTATIONS_BASE_NAME
                phi_file = os.path.join(actual_angles_root, phi_name)
                angles_name = ANGLES_BASE_NAME + str(orientation) + '.mat'
                angles_file = os.path.join(actual_angles_root, angles_name)
                chi_name = CHI_BASE_NAME + str(orientation) + ORIENTATIONS_BASE_NAME
                mms_name = MMS_BASE_NAME + str(orientation) + ORIENTATIONS_BASE_NAME
                msa_name = MSA_BASE_NAME + str(orientation) + ORIENTATIONS_BASE_NAME
                v1_name = V1_BASE_NAME + str(orientation) + ORIENTATIONS_BASE_NAME
                chi_reconstructed_file = os.path.join(actual_angles_root, chi_name)
                mms_file = os.path.join(actual_angles_root, mms_name)
                msa_file = os.path.join(actual_angles_root, msa_name)
                v1_file = os.path.join(actual_angles_root, v1_name)
                mms_diff_name = os.path.join(actual_angles_root, 'diff_mms.nii.gz')
                msa_diff_name = os.path.join(actual_angles_root, 'diff_msa.nii.gz')
                v1_diff_name = os.path.join(actual_angles_root, 'diff_v1.nii.gz')

                check_root(SIMULATED_FOLDER)
                check_root(actual_theta_folder)
                check_root(actual_angles_root)

                print('Reading sti image')
                chi_name = CHI_ROOT + CHI_NAMES
                chi, mat_size, fov, nii = read_img(chi_name)
                mask, _, _, _, = read_img(MASK_FILE)
                print('Done')
                print('---------------------------------------')
                print('Simulating' + str(orientation) + ' orientations')
                print('... Calculating off-resonant phase')
                vec_theta, vec_psi = angles_cylinders(orientation, n_theta, n_psi)
                mat_projection = projection_variables(vec_theta, vec_psi, orientation, fov, mat_size)
                phi = get_total_phase(chi, mat_projection)
                phi = phase_scale * (phi + STD_NOISE * torch.rand(phi.size()))
                phi = mask_image(phi, mask)
                print('Saving phase to nifti')
                save_img(nii, phi, phi_file)
                dic = {'vec_theta': vec_theta.numpy(), 'vec_psi': vec_psi.numpy()}
                savemat(angles_file, dic)
                chi_rec = min_squares(phi / phase_scale, mat_projection)
                print(chi_rec.shape)
                print(mask.shape)
                chi_rec = mask_image(chi_rec, mask)
                print('... Saving reconstructed tensor to nifti')
                save_img(nii, chi_rec, chi_reconstructed_file)

                print('Generating eigen-decomposition and calculating susceptibility components (mms and msa)')
                v1, mms, msa = gen_mms_msa(chi_rec)
                print('... Saving eigen components')
                save_img(nii, mms, mms_file)
                save_img(nii, msa, msa_file)
                save_img(nii, v1, v1_file)

                tmp = mms - mms_gt
                save_img(nii, tmp, mms_diff_name)
                tmp = msa - msa_gt
                save_img(nii, tmp, msa_diff_name)
                tmp = get_angle_images(v1_gt, v1)
                save_img(nii, tmp * white_matter, v1_diff_name)
                print('Done')


if __name__ == '__main__':
    main()
