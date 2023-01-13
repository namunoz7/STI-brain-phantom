import tensorflow as tf
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import os

ROOT_DIR = os.path.join('..', '..', 'Imagenes', 'Phantom_real_data')
ROOT_MODEL = os.path.join(ROOT_DIR, 'Phantom_tensor')
RECONSTRUCTION_ROOT = os.path.join(ROOT_DIR, 'Simulated_data')
ROOT_DIFFUSION = os.path.join(RECONSTRUCTION_ROOT, 'Diffusion_sti')
ROOT_SUSCEPTIBILITY = os.path.join(RECONSTRUCTION_ROOT, 'Susceptibility_sti')
ROOT_6_ORIENTATIONS = '6_orientations'
IMG_ROOT = os.path.join(RECONSTRUCTION_ROOT, 'Phantom_no_microstructure')
CHI_DIFFUSION_NAME = os.path.join(ROOT_MODEL, 'chi_dti_filt.nii.gz')
CHI_SUSCEPTIBILITY_NAME = os.path.join(ROOT_MODEL, 'chi_sti_filt.nii.gz')
CHI_1_NAME = os.path.join(IMG_ROOT, 'chi_6_orientations.nii.gz')
CHI_2_NAME = os.path.join(IMG_ROOT, 'chi_12_orientations.nii.gz')
MASK_NAME = os.path.join(ROOT_DIR, 'Masks', 'brain_mask.nii.gz')

phantom_models = [CHI_DIFFUSION_NAME, CHI_SUSCEPTIBILITY_NAME]


def read_img(filename):
    """
    Read nifty image
    :param filename:
    :return:
    """
    nii_img = nib.load(filename)
    img = nii_img.get_fdata()
    # img_shape = np.array(img.shape[0:3])
    # header = nii_img.header
    # voxel_size = np.array(header.get_zooms()[0:3])
    # fov = img_shape * voxel_size

    return img  # , img_shape, fov, nii_img


def image_tensor(tensor_1, tensor_2, img, title_name):
    """
    Returns a 3x3 grid susceptibility images, displaying the susceptibility tensor as a matplotlib figure
    :param tensor_2: Susceptibility tensor from susceptibility eigenvector
    :param title_name: Title name of the figure
    :param img: Slice of the images to plot
    :param tensor_1: Susceptibility tensor from diffusion eigenvector
    :return: figure
    """
    # Create a figure containing the plot
    fig, axs = plt.subplots(2, 3, figsize=(11.3, 8))
    indexes = [0, 3, 5]
    positions = [0.03, 0.33, 0.63, 0.03, 0.33, 0.63]
    fig.suptitle(title_name, fontsize=20)
    legends = [r'$\chi_{11}$', r'$\chi_{22}$', r'$\chi_{33}$']
    tmp_pos = 0
    for n in range(2):
        if n == 0:
            tensor = tensor_1
        else:
            tensor = tensor_2
        for m in range(3):
            if n == 0:
                tmp = [positions[tmp_pos], 0.53, 0.35, 0.43]
            else:
                tmp = [positions[tmp_pos], 0.1, 0.35, 0.43]
            # Start next subplot
            axs[n, m].imshow(np.fliplr(np.rot90(tensor[:, :, img, indexes[m]])), cmap='gray', vmin=-0.1, vmax=0.1)
            axs[n, m].set_position(pos=tmp)
            axs[n, m].text(20, 40, legends[m], color='white', fontsize=18)
            axs[n, m].set_xticks([])
            axs[n, m].set_yticks([])
            tmp_pos += 1
    mappable = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=-0.1, vmax=0.1), cmap='gray')
    cax = plt.axes([0.04, 0.05, 0.93, 0.03])
    bar = plt.colorbar(mappable=mappable, cax=cax, orientation='horizontal')
    bar.ax.tick_params(labelsize=13)
    return fig


def imshow_sti(tensor, img, title_name):
    """
    Returns a 3x3 grid susceptibility images, displaying the susceptibility tensor as a matplotlib figure
    :param title_name: Title name of the figure
    :param img: Slice of the images to plot
    :param tensor: Susceptibility image tensor to be displayed
    :return: figure
    """
    # Create a figure containing the plot
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 9.5))
    indexes = [0, 3, 5, 1, 2, 4]
    positions = [0.03, 0.33, 0.63, 0.03, 0.33, 0.63]
    fig.suptitle(title_name, fontsize=20)
    legends = [r'$\chi_{11}$', r'$\chi_{22}$', r'$\chi_{33}$', r'$\chi_{12}$', r'$\chi_{13}$', r'$\chi_{23}$']
    tmp_pos = 0
    for n in range(2):
        for m in range(3):
            if n == 0:
                tmp = [positions[tmp_pos], 0.53, 0.35, 0.43]
            else:
                tmp = [positions[tmp_pos], 0.1, 0.35, 0.43]
            # Start next subplot
            axs[n, m].imshow(np.fliplr(np.rot90(tensor[:, img, :, indexes[tmp_pos]])), cmap='gray', vmin=-0.1, vmax=0.1)
            axs[n, m].set_position(pos=tmp)
            axs[n, m].text(20, 40, legends[tmp_pos], color='white', fontsize=18)
            axs[n, m].set_xticks([])
            axs[n, m].set_yticks([])
            tmp_pos += 1
    mappable = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=-0.1, vmax=0.1), cmap='gray')
    cax = plt.axes([0.04, 0.05, 0.93, 0.03])
    bar = plt.colorbar(mappable=mappable, cax=cax, orientation='horizontal')
    bar.ax.tick_params(labelsize=13)
    return fig


def image_phi(tensor, img):
    """
    Returns a 3x3 grid susceptibility images, displaying the susceptibility tensor as a matplotlib figure
    :param img: Slice number of the image
    :param tensor: Susceptibility image tensor to be displayed
    :return: figure
    """
    # Create a figure containing the plot
    figure = plt.figure()
    for n in range(0, 8):
        # Start next subplot
        plt.subplot(3, 3, n+1)
        plt.xticks([])
        plt.yticks([])
        plt.imshow(np.fliplr(np.rot90(tensor[:, :, img, n])), cmap='gray', vmin=-0.05, vmax=0.05)
    cax = plt.axes([0.12, 0.05, 0.78, 0.05])
    plt.colorbar(cax=cax, orientation='horizontal')
    return figure


def hist_tensor(tensor):
    """
    Plots histogram of the input tensor
    :param tensor:
    :return:
    """
    iso_chi = tf.concat([tensor[:, 0], tensor[:, 3], tensor[:, 5]], 0)
    ani_chi = tf.concat([tensor[:, 1], tensor[:, 2], tensor[:, 4]], 0)
    # img = tensor[:, n_img]
    # img = tf.reshape(img, (1, tf.math.reduce_prod(tf.size(img))))
    fig = plt.figure()
    plt.hist(iso_chi.numpy(), range=(-0.1, 0.1), bins=100, label='Isotropic')
    plt.hist(ani_chi.numpy(), range=(-0.1, 0.1), bins=100, label='Anisotropic', alpha=0.7)
    plt.legend()
    fig.suptitle('Phantom Distribution', fontsize=18)
    return fig


def hist_eig_values(tensor, mask):
    """
    Plots histogram of the input tensor
    :param mask:
    :param tensor:
    :return:
    """
    tensor = tf.boolean_mask(tensor, mask)
    eig1, eig2, eig3 = tf.split(tensor, num_or_size_splits=3, axis=1)
    fig = plt.figure()
    plt.hist(tf.reshape(eig1, shape=tf.size(eig1)).numpy(),
             range=(-0.1, 0.1), bins=200, label='Eigenvalue 1')
    plt.hist(tf.reshape(eig2, shape=tf.size(eig2)).numpy(),
             range=(-0.1, 0.1), bins=200, label='Eigenvalue 2', alpha=0.85)
    plt.hist(tf.reshape(eig3, shape=tf.size(eig3)).numpy(),
             range=(-0.1, 0.1), bins=200, label='Eigenvalue 3', alpha=0.7)
    plt.legend()
    fig.suptitle('Phantom eigenvalues', fontsize=18)
    return fig


def main():
    chi_dti = read_img(CHI_DIFFUSION_NAME)
    chi_sti = read_img(CHI_SUSCEPTIBILITY_NAME)
    image_tensor(chi_dti, chi_sti, 46, '')
    plt.show()
    print('Done')


if __name__ == '__main__':
    main()
