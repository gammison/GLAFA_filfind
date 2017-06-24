'''
tree visualization lib
'''

import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

fil_finder_dir = '/Users/larryli/Documents/CC/16-17/research/GALFA_filfind'
sys.path.append(fil_finder_dir)

import mask_obj_node_tree as maskTree
import mask_obj_node as maskNode

import cube_processes as cube


def vis_tree_shadow(tree, tree_name, save_fig=False, save_dir=None, save_name=None):
    root = tree.root_node
    starting_v = tree.root_v_slice
    length = tree.length
    assert(tree.has_ended)
    plot_corners = [root.corner_BL[0], root.corner_TR[0], root.corner_BL[1], root.corner_TR[1]]

    if save_fig and save_name is None:
        save_name = tree_name

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow([[], []], extent=plot_corners)
    ax.contour(root.mask, colors='r', extent=plot_corners)

    fig.suptitle(tree_name)
    ax.set_xlabel('tree length: {}'.format(length))

    if save_fig:
        fig.savefig(save_dir + save_name)
