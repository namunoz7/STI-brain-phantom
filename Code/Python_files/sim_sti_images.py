import scipy.io
import torch
import nibabel as nib
import numpy as np
import os
import sys
import scipy
from datetime import datetime

# IMG_ROOT = os.path.join('..', '..', '..', 'Imagenes', 'Phantom_real_data')
IMG_ROOT = os.path.join('..', 'Phantom_real_data')
MASK_FILES = os.path.join(IMG_ROOT, 'Masks', 'brain_mask.nii.gz')
MODEL_TENSOR = ['Diffusion_sti', 'Susceptibility_sti']
TMP_ORIENTATIONS = ['Simulating 6 orientations', 'Simulating 12 orientations']
DIR_ORIENTATIONS = ['6_orientations', '12_orientations']
PHASE_FILES = ['phi_6_orientations.nii.gz', 'phi_12_orientations.nii.gz']
CHI_FILE = 'chi.nii.gz'
ANGLES_FILES = ['angles_1.mat', 'angles_2.mat']
N_ORIENTATIONS = [6, 12]
EIG_VAL_NAME = 'eig_val.nii.gz'
EIG_V1_NAME = 'eig_v1.nii.gz'
EIG_V2_NAME = 'eig_v2.nii.gz'
EIG_V3_NAME = 'eig_v3.nii.gz'
MMS_NAME = 'mms.nii.gz'
MSA_NAME = 'msa.nii.gz'

eps = sys.float_info.epsilon
GAMMA = 42.58
TE = 1e-3
B0 = 7  # T
phase_scale = 2*np.pi*GAMMA*B0*TE


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

    k = torch.stack([kxx, kyy, kzz], dim=-1)
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

    fft_phi = fft_phase(bulk_phase)
    tmp_real = torch.matmul(mat_transpose, torch.real(fft_phi).unsqueeze(-1))
    tmp_img = torch.matmul(mat_transpose, torch.imag(fft_phi).unsqueeze(-1))

    real_ft_chi = torch.matmul(b_inv, tmp_real).unsqueeze(-1)
    img_ft_chi = torch.matmul(b_inv, tmp_img).unsqueeze(-1)
    tmp_chi = torch.cat((real_ft_chi, img_ft_chi), dim=-1)
    chi = inv_fft_phase(torch.view_as_complex(tmp_chi))
    chi = torch.real(chi)
    return torch.squeeze(chi)


def save_img(header, img, name):
    """
    Save the image as a nifti file
    :param header:
    :param img:
    :param name:
    """
    new_dim = (img.size(0), img.size(1), img.size(2), img.size(-1))
    header.set_data_shape(new_dim)
    nii_phi = nib.Nifti1Image(img, affine=None, header=header)
    nib.save(nii_phi, name)
    print('...... ' + name + ' saved')


def sim_tensor(phase_file, angles_file, num_orientations, phi_scale, mask):
    """
    Simulates the susceptibility tensor, according to the local phase of each orientation, the number of orientations
    simulated, and saves the tensor in nifti format.
    :param phi_scale: Phase scale to apply to the phase.
    :param mask: Boolean mask of teh tensor
    :param phase_file: Name of the phase simulated. It must be saved in nifti format (.nii.gz)
    :param angles_file: Angles that were simulated. Must correspond to an array of 2 columns, and number of rows
    depending on the number of orientations simulated.
    :param num_orientations: Number of orientations simulated.
    :return:
    """
    print('...... Reading images')
    phi, mat_size, fov, nii = read_img(phase_file)
    mat = scipy.io.loadmat(angles_file)
    vec_theta = torch.tensor(mat['vec_theta'])
    vec_psi = torch.tensor(mat['vec_psi'])
    mask, _, _, _ = read_img(mask)
    mask = mask.unsqueeze(-1).repeat([1, 1, 1, 6])
    print('...... Calculating Susceptibility Tensor')
    chi = min_squares(phi/phi_scale, vec_theta, vec_psi, num_orientations, fov, mat_size)
    chi = torch.mul(chi, mask)
    return chi


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
    l, v = torch.linalg.eigh(tensor)
    v1 = v[:, :, :, :, 2]
    v2 = v[:, :, :, :, 1]
    v3 = v[:, :, :, :, 0]
    return l, v1, v2, v3


def gen_mms_msa(eigenvalues):
    """
    Calculates the mean magnetic susceptibility (MMS) and magnetic susceptibility anisotropy (MSA)
    :param eigenvalues: Susceptibility tensor input in vectorized form
    :return: mms, msa
    """
    mms = torch.mean(eigenvalues, dim=-1)
    msa = eigenvalues[:, :, :, 2] - torch.mean(eigenvalues[:, :, :, 0:2], dim=-1)
    return mms, msa


def main():
    for n_mode in MODEL_TENSOR:
        actual_model = os.path.join(IMG_ROOT, 'Simulated_data', n_mode)
        print('====================================')
        print(n_mode)
        for n_orientation in range(2):
            print(TMP_ORIENTATIONS[n_orientation])
            actual_folder = os.path.join(actual_model, DIR_ORIENTATIONS[n_orientation])
            phi_file = os.path.join(actual_folder, PHASE_FILES[n_orientation])
            angles_file = os.path.join(actual_folder, ANGLES_FILES[n_orientation])
            sti_name = os.path.join(actual_model, CHI_FILE)
            eig_name = os.path.join(actual_model, EIG_VAL_NAME)
            v1_name = os.path.join(actual_model, EIG_V1_NAME)
            v2_name = os.path.join(actual_model, EIG_V2_NAME)
            v3_name = os.path.join(actual_model, EIG_V3_NAME)
            mms_name = os.path.join(actual_folder, MMS_NAME)
            msa_name = os.path.join(actual_folder, MSA_NAME)

            _, _, _, nii = read_img(phi_file)
            chi = sim_tensor(phase_file=phi_file, angles_file=angles_file,
                             num_orientations=N_ORIENTATIONS[n_orientation], phi_scale=phase_scale, mask=MASK_FILES)
            print('... Saving tensor')
            hdr = nii.header
            new_dim = (chi.size(0), chi.size(1), chi.size(2), chi.size(-1))
            hdr.set_data_shape(new_dim)
            save_img(hdr, chi, sti_name)
            print('...  Eigen-decomposition')
            l, v1, v2, v3 = get_eigen_decomposition(chi)
            print('... Saving eigen-decomposition')
            save_img(hdr, l, eig_name)
            save_img(hdr, v1, v1_name)
            save_img(hdr, v2, v2_name)
            save_img(hdr, v3, v3_name)
            mms, msa = gen_mms_msa(l)
            print('... Saving MMS and MSA')
            new_dim = (mms.size(0), mms.size(1), mms.size(2))
            hdr.set_data_shape(new_dim)
            nii = nib.Nifti1Image(mms, affine=None, header=hdr)
            nib.save(nii, mms_name)
            print('... ' + mms_name + ' saved.')
            nii = nib.Nifti1Image(msa, affine=None, header=hdr)
            nib.save(nii, msa_name)
            print('... ' + msa_name + ' saved.')
            print('----------------------------')

    print(datetime.now())


if __name__ == '__main__':
    main()
