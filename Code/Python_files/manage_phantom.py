import torch
import nibabel as nib
import numpy as np
import sys
import os
from scipy.io import savemat

# IMG_ROOT = '../../Imagenes/Phantom_real_data/'
IMG_ROOT = '../Phantom_real_data/'
# IMG_ROOT = 'Images/Data_Phantom/data/stimodel/'
# MASK_FILE = 'Images/Data_Phantom/data/masks/mask_tensor3.nii.gz'

MASK_FILE = IMG_ROOT + 'Masks/brain_mask.nii.gz'
SIMULATED_FOLDER = IMG_ROOT + 'Noised_data/'
# SIMULATED_FOLDER = IMG_ROOT + 'Simulated_data/'
ORIENTATION_6_ROOT = '6_orientations/'
ORIENTATION_12_ROOT = '12_orientations/'
CHI_ROOT = IMG_ROOT + 'Phantom_tensor/'
CHI_NAMES = ['chi_dti_filt.nii.gz', 'chi_sti_filt.nii.gz']
SIM_TENSOR_FOLDER = ['Diffusion_sti/', 'Susceptibility_sti/']

PHI_NAME_1 = 'phi_6_orientations.nii.gz'
CHI_RECONSTRUCTED_1 = 'chi_6_orientations.nii.gz'
ANGLES_NAME_1 = 'angles_1.mat'
N_ORIENTATIONS_1 = 6

PHI_NAME_2 = 'phi_12_orientations.nii.gz'
CHI_RECONSTRUCTED_2 = 'chi_12_orientations.nii.gz'
ANGLES_NAME_2 = 'angles_2.mat'
N_ORIENTATIONS_2 = 12

eps = sys.float_info.epsilon
GAMMA = 42.58
TE = 12e-3
B0 = 3  # T
# phase_scale = 2*np.pi*GAMMA*TE*B0
phase_scale = 2*np.pi*GAMMA*B0

max_theta_1 = 50.0
max_phi_1 = 50.0
max_theta_2 = 50.0
max_phi_2 = 50.0
std_noise = 1e-6


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


def angles_cylinders(num_rotations, max_theta, max_phi):
    """
    Rotation angles that are used to rotate the object in the scanner
    :param num_rotations:
    :param max_theta:
    :param max_phi:
    """
    # Parameters of the angles to model the cylinders
    theta_min = 10.0
    theta2 = max_theta
    psi2 = max_phi

    # Tilt angle of main field
    vec_theta = (theta2 - theta_min) * torch.rand(num_rotations-1, dtype=torch.float64) + theta_min
    vec_theta = vec_theta * np.pi / 180.0
    vec_theta = torch.cat((torch.tensor((0.0,), dtype=torch.float64), vec_theta), 0)
    vec_theta = vec_theta.reshape(1, num_rotations)

    # Rotation angle with the z axis
    vec_psi = 2 * psi2 * torch.rand(num_rotations-1, dtype=torch.float64) - psi2
    vec_psi = vec_psi * np.pi / 180.0
    vec_psi = torch.cat((torch.tensor((0.0,), dtype=torch.float64), vec_psi), 0)
    vec_psi = vec_psi.reshape(1, num_rotations)

    return vec_theta, vec_psi


def get_direction_field(vec_theta, vec_psi, n_orientations):
    """
    Gets the direction field vector of the multiple orientations, made by the cylinders. All the angles are in radians
    :param vec_theta: tilt angle (inclination angle with respect to the z axis)
    :param vec_psi: azimuth angle (rotation angle made in the x-y plane)
    :param n_orientations:
    :return:
    """
    direction_field = torch.zeros(n_orientations, 3, dtype=torch.float64)
    direction_field[:, 0] = torch.sin(vec_theta) * torch.cos(vec_psi)  # Hx
    direction_field[:, 1] = torch.sin(vec_theta) * torch.sin(vec_psi)  # Hy
    direction_field[:, 2] = torch.cos(vec_theta)  # Hz
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


def get_total_phase(chi, vec_theta, vec_psi, fov, mat_size, n_orientations):
    """
    Calculates the total phase of the tensor using the linear form of the Tensor model
    :param n_orientations: Number of orientations in the
    :param mat_size:
    :param fov:
    :param chi: Vectorized form of the tensor
    :param vec_theta: Angle of deviation of the image with the main magnetic field direction
    :param vec_psi: Angle of rotation in the x-y plane
    :return: phase: the local phase of the image
    """
    print('... Generating matrix projection')
    matrix_projection = projection_variables(vec_theta, vec_psi, n_orientations, fov, mat_size)
    print('... Done')
    print('... Calculating off-resonant field')
    print('...... FFT of susceptibility tensor')
    fft_chi = fft_phase(chi)
    print('...... Multiplication of matrix projection and FFT of susceptibility tensor')
    tmp_real = torch.matmul(matrix_projection, torch.real(fft_chi).unsqueeze(-1))
    tmp_img = torch.matmul(matrix_projection, torch.imag(fft_chi).unsqueeze(-1))
    print('...... Inverse FFT of off-resonant field')
    tmp_phi = torch.cat((tmp_real, tmp_img), dim=-1)
    phi = inv_fft_phase(torch.view_as_complex(tmp_phi))
    return torch.real(phi), matrix_projection


def min_squares(bulk_phase, vec_theta, vec_psi, n_orientations, fov, mat_size):
    """
    Calculates the tensor by getting the phase of geometric figures.
    :param vec_theta:
    :param vec_psi:
    :param n_orientations:
    :param fov:
    :param mat_size:
    :param bulk_phase: Total phase of the image phase
    :return: The vectorized form of the tensor
    """
    # import Generate_Dataset.fourier_torch as ft

    print('... Generating matrix projection')
    vec_theta = torch.round(vec_theta*180.0/np.pi)*np.pi/180.0
    vec_psi = torch.round(vec_psi * 180.0 / np.pi) * np.pi / 180.0
    mat_projection = projection_variables(vec_theta, vec_psi, n_orientations, fov, mat_size)
    mat_transpose = mat_projection.transpose(3, 4)
    b_mat = torch.matmul(mat_transpose, mat_projection)
    b_inv = b_mat.inverse()
    print('b_inv shape: ' + str(b_inv.shape))

    fft_phi = fft_phase(bulk_phase)
    print('FFT phi shape: ' + str(fft_phi.shape))
    tmp_real = torch.matmul(mat_transpose, torch.real(fft_phi).unsqueeze(-1))
    print('FFT tmp_real shape: ' + str(tmp_real.shape))
    tmp_img = torch.matmul(mat_transpose, torch.imag(fft_phi).unsqueeze(-1))
    print('FFT tmp_img shape: ' + str(tmp_img.shape))
    real_ft_chi = torch.matmul(b_inv, tmp_real).unsqueeze(-1)
    print('FFT real chi: ' + str(real_ft_chi.shape))
    img_ft_chi = torch.matmul(b_inv, tmp_img).unsqueeze(-1)
    print('FFT imaginary chi: ' + str(img_ft_chi.shape))

    tmp_chi = torch.cat((real_ft_chi, img_ft_chi), dim=-1)
    print('FFT chi: ' + str(tmp_chi.shape))
    chi = inv_fft_phase(torch.view_as_complex(tmp_chi))
    return torch.real(chi)


def save_img(nii, img, name):
    """
    Save the image as a nifti file
    :param nii:
    :param img:
    :param name:
    """
    header = nii.header
    nii_phi = nib.Nifti1Image(img.numpy(), header.get_best_affine())
    nib.save(nii_phi, name)


def main():
    vec_theta_6, vec_psi_6 = angles_cylinders(N_ORIENTATIONS_1, max_theta_1, max_phi_1)
    vec_theta_12, vec_psi_12 = angles_cylinders(N_ORIENTATIONS_2, max_theta_2, max_phi_2)
    for n_chi in range(2):
        print('=======================================')
        actual_model_root = SIMULATED_FOLDER + SIM_TENSOR_FOLDER[n_chi]
        if not os.path.exists(actual_model_root):
            os.mkdir(actual_model_root)
            print("Directory ", actual_model_root, " Created ")
        else:
            print("Directory ", actual_model_root, " already exists")

        if n_chi == 0:
            print('Using STI with diffusion eigenvectors')
        else:
            print('Using STI with susceptibility eigenvectors')
        print('Reading sti image')
        chi_name = CHI_ROOT + CHI_NAMES[n_chi]
        actual_orientation = actual_model_root + ORIENTATION_6_ROOT
        phi_file = actual_orientation + PHI_NAME_1
        angles_file = actual_orientation + ANGLES_NAME_1
        chi_reconstructed_file = actual_orientation + CHI_RECONSTRUCTED_1

        chi, mat_size, fov, nii = read_img(chi_name)
        mask, _, _, _, = read_img(MASK_FILE)
        print('Done')

        print('---------------------------------------')
        if not os.path.exists(actual_orientation):
            os.mkdir(actual_orientation)
            print("Directory ", actual_orientation, " Created ")
        else:
            print("Directory ", actual_orientation, " already exists")
        print('vec_theta_1:')
        print(vec_theta_6)
        print('vec psi_1:')
        print(vec_psi_6)
        print('Calculating off-resonant phase')
        phi, mat_projection = get_total_phase(chi, vec_theta_6, vec_psi_6, fov, mat_size, N_ORIENTATIONS_1)
        phi = torch.mul(phi*phase_scale, mask.unsqueeze(-1).repeat([1, 1, 1, N_ORIENTATIONS_1]))
        phi = phi/phase_scale + std_noise * torch.rand(phi.size())
        print('Done')
        print('Saving phase to nifti')
        save_img(nii, phi, phi_file)
        print('Done')
        dic = {'vec_theta': vec_theta_6.numpy(), 'vec_psi': vec_psi_6.numpy()}
        savemat(angles_file, dic)
        chi_1 = min_squares(phi, vec_theta_6, vec_psi_6, N_ORIENTATIONS_1, fov, mat_size)
        print('Saving reconstructed tensor to nifti')
        save_img(nii, chi_1, chi_reconstructed_file)

        print('---------------------------------------')
        actual_orientation = actual_model_root + ORIENTATION_12_ROOT
        phi_file = actual_orientation + PHI_NAME_2
        angles_file = actual_orientation + ANGLES_NAME_2
        chi_reconstructed_file = actual_orientation + CHI_RECONSTRUCTED_2
        if not os.path.exists(actual_orientation):
            os.mkdir(actual_orientation)
            print("Directory ", actual_orientation, " Created ")
        else:
            print("Directory ", actual_orientation, " already exists")
        print('Simulating 6 orientations')
        print('Simulating 12 orientations')
        print('vec theta 2:')
        print(vec_theta_12)
        print('vec psi 2:')
        print(vec_psi_12)
        print('Calculating off-resonant phase')
        phi, mat_projection = get_total_phase(chi, vec_theta_12, vec_psi_12, fov, mat_size, N_ORIENTATIONS_2)
        phi = torch.mul(phi*phase_scale, mask.unsqueeze(-1).repeat([1, 1, 1, N_ORIENTATIONS_2]))
        phi = phi / phase_scale + std_noise * torch.rand(phi.size())
        print('Done')
        print('Saving phase to nifti')
        save_img(nii, phi, phi_file)
        print('Done')
        dic = {'vec_theta': vec_theta_12.numpy(), 'vec_psi': vec_psi_12.numpy()}
        savemat(angles_file, dic)
        chi_2 = min_squares(phi, vec_theta_12, vec_psi_12, N_ORIENTATIONS_2, fov, mat_size)
        print('Saving reconstructed tensor to nifti')
        save_img(nii, chi_2, chi_reconstructed_file)


if __name__ == '__main__':
    main()
