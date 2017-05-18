'''
contains all the functions related to processing a cube and parsing data

'''

import sys
import pickle
import scipy.ndimage
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as coord
import mask_obj_node_tree as maskTree
import mask_obj_node as maskNode


def find_all_trees_from_slices(vs, dict_full_names, hdr, overlap_thresh=.75, reverse_find=False, verbose=False):
    '''
    find all the trees given a list of velocity and a list of dict names
    reverse find currently not implemented!

    Arguments:
        vs {[type]} -- [description]

    Keyword Arguments:
        overlap_thresh {number} -- [description] (default: {.75})
        reverse_find {bool} -- [description] (default: {False})
        verbose {bool} -- [description] (default: {False})
    '''
    nodes_by_tree = {}

    if len(vs) != len(dict_full_names):
        print "\n\nLength of velocities and dictionaries don't match"
        sys.exit()

    for v in range(len(vs)):
        dict_recover_path = dict_full_names[v]
        nodes_in_v_slice = pickle.load(open(dict_recover_path, 'rb'))

        if verbose:
            print "working on v slice %d" % vs[v]

        if v == 0:
            for j in sorted(nodes_in_v_slice.keys(), reverse=True):
                # create trees and match masks by descending size
                if verbose:
                    print "\ton mask %d" % j

                current_node = nodes_in_v_slice[j]
                if not maskNode.check_node_b_cutoff(current_node, hdr):
                    continue

                new_tree = maskTree.newTreeFromNode(current_node, verbose=verbose)
                add_tree_to_dict(new_tree, nodes_by_tree)
        else:
            for j in sorted(nodes_in_v_slice.keys(), reverse=True):
                # match masks by descending size onto existing trees
                if verbose:
                    print "\ton mask %d" % j

                current_node = nodes_in_v_slice[j]
                if not maskNode.check_node_b_cutoff(current_node, hdr):
                    continue

                if not match_node_onto_tree(current_node, nodes_by_tree, overlap_thresh, verbose=verbose):
                    new_tree = maskTree.newTreeFromNode(current_node, verbose=verbose)
                    add_tree_to_dict(new_tree, nodes_by_tree)

            end_noncontinuous_trees(nodes_by_tree, vs[v])
            delete_small_dead_trees(nodes_by_tree)

        del nodes_in_v_slice

    return nodes_by_tree


def match_node_onto_tree(node, trees, overlap_thresh, verbose=False):
    '''
    matches a node onto an existing tree

    Arguments:
        node {[type]} -- [description]
        tree {[type]} -- [description]
    '''
    has_matched = False

    for k in sorted(trees.keys(), reverse=True):
        # search the existing trees by descending size
        if trees[k].has_ended:
            continue
        else:
            if trees[k].getLastNode().checkMaskOverlap(node, overlap_thresh):
                node.visited = True
                has_matched = True
                trees[k].addNode(node, verbose)
                if verbose:
                    print "\t\t matched with tree %s" % k
                break
            else:
                continue

    return has_matched


def end_noncontinuous_trees(trees, current_v):
    '''
    ends the trees that have not matched in this velocity slice

    Arguments:
        trees {[type]} -- [description]
        current_v {[type]} -- [description]
    '''
    for i in sorted(trees.keys()):
        if trees[i].has_ended:
            continue
        else:
            if trees[i].root_v_slice + trees[i].length <= current_v:
                trees[i].has_ended = True


def find_all_trees_from_cube(nodes_by_v_slice, overlap_thresh=.75, reverse_find=False, verbose=False):
    '''
    find all the trees given the entire set of dicts which contain the entire
    cube
    '''
    nodes_by_tree = {}

    vs = sorted(nodes_by_v_slice.keys(), reverse=reverse_find)

    for v in vs:
        if v not in nodes_by_v_slice:
            print "\n\nSOMETHING WENT WRONG"
            sys.exit()

        if verbose:
            print "working on v slice %d" % v

        for i in sorted(nodes_by_v_slice[v].keys(), reverse=True):
            # create trees and match masks by descending size
            if verbose:
                print "\ton mask %d" % i

            current_node = nodes_by_v_slice[v][i]

            if current_node.visited == True:
                if verbose:
                    print "\tvisited"
                continue
            else:
                new_full_tree = find_new_full_tree(current_node, nodes_by_v_slice, overlap_thresh, reverse=reverse_find)
                nodes_by_tree[new_full_tree.getTreeMaskedArea2D()] = new_full_tree

    return nodes_by_tree


def find_new_full_tree(root_node, nodes_by_v_slice, overlap_thresh, reverse=False, verbose=False):
    new_tree = maskTree.newTreeFromNode(root_node, verbose=verbose)
    new_tree = find_all_children(new_tree, nodes_by_v_slice, overlap_thresh, reverse=reverse, verbose=verbose)

    if verbose:
        print new_tree.root_node.corners

    return new_tree


def find_all_children(tree, nodes_by_v_slice, overlap_thresh, reverse=False, verbose=False):
    start_v = tree.root_v_slice

    if reverse:
        adv = -1
    else:
        adv = 1

    v_slice = start_v + adv
    while True:
        if v_slice not in nodes_by_v_slice or len(nodes_by_v_slice[v_slice]) == 0:
            tree.has_ended = True
            break
        else:
            children_added = []
            if verbose:
                print "Corners of old nodes" + str(tree.getLastNode().corners)

            for i in sorted(nodes_by_v_slice[v_slice].keys(), reverse=True):
                # matches masks by big to small
                this_node = nodes_by_v_slice[v_slice][i]
                if this_node.visited == False:
                    if tree.getLastNode().checkMaskOverlap(this_node, overlap_thresh):
                        children_added.append(i)
                        if verbose:
                            print "Corners of matched children" + str(nodes_by_v_slice[v_slice][i].corners)

            if verbose:
                print "\t\t %d children found on slice %d" % (len(children_added), v_slice)

            if len(children_added) == 1:
                nodes_by_v_slice[v_slice][children_added[0]].visited = True
                if verbose:
                    print "\t\t\t %d - %d marked as visited" % (v_slice, children_added[0])
                tree.addNode(nodes_by_v_slice[v_slice][children_added[0]])
                v_slice += adv
                continue
            elif len(children_added) > 1:
                first_node = nodes_by_v_slice[v_slice][children_added[0]]
                first_node.visited = True
                for j in children_added:
                    if j != children_added[0]:
                        first_node.mergeNode(nodes_by_v_slice[v_slice][j])
                        nodes_by_v_slice[v_slice][j].visited = True
                    if verbose:
                        print "\t\t\t %d - %d marked as visited" % (v_slice, children_added[j])
                tree.addNode(first_node)
                v_slice += adv
                continue
            elif len(children_added) == 0:
                tree.has_ended = True
                break

    return tree


def prune_trees(all_trees, size_cut=30, length_cut=3, length_limit=36, verbose=False):
    pruned_trees = dict(all_trees)

    bad_trees_keys = []

    for k in pruned_trees:
        # cut length
        if pruned_trees[k].length >= length_limit or pruned_trees[k].length < length_cut:
            bad_trees_keys.append(k)
            continue

        # cut size
        if pruned_trees[k].getTreeMaskedArea2D() < size_cut:
            bad_trees_keys.append(k)
            continue

    if verbose:
        print "\t\t deleting %d trees that don't fit the criteria" % len(bad_trees_keys)

    for k in bad_trees_keys:
        del pruned_trees[k]

    return pruned_trees


def get_length_dist(all_trees):
    '''
    return the length distribution of all the trees
    Arguments:
        all_trees {[type]} -- [description]

    Returns:
        [type] -- [description]
    '''
    length_dist = []

    for k in all_trees:
        length_dist.append(all_trees[k].length)

    return np.array(length_dist)


def index_to_radec(xs, ys, hdr, verbose=True):
    '''
    turns arrays of indecies to ra, dec based on header

    Arguments:
        xs {[type]} -- [description]
        ys {[type]} -- [description]
        hdr {[type]} -- [description]
    '''

    if type(xs) != np.ndarray:
        xs = np.array(xs)
    if type(ys) != np.ndarray:
        ys = np.array(ys)

    if xs.size != ys.size:
        if verbose:
            print "\t x & y array sizes different in conversion to RA-DEC"

    ras = (xs - hdr['CRPIX1']) * hdr['CDELT1'] + hdr['CRVAL1']
    decs = (ys - hdr['CRPIX2']) * hdr['CDELT2'] + hdr['CRVAL2']
    # ^redundant but the right way?
    return ras, decs


def mask_lb(data, hdr, b_cutoff=30, toNaN=False, l_cutoff=None):
    pass


def umask_and_save(data, hdr, save_dir, file_name, radius=None):
    '''
    unsharp masking a data slice and saving it
    '''
    save_path = save_dir + file_name.rsplit('.', 1)[0] + '_umask.fits'
    if radius is None:
        umask_data = umask(data)
    else:
        umask_data = umask(data, radius)
    fits.writeto(save_path, umask_data, header=hdr)

    return umask_data


def radecs_to_lb(ras, decs):
    '''
    Transformation between lists of ras, decs, to ls, bs. Assumes ra, dec in degrees
    Conforms to astropy 0.4.3
    taken from https://github.com/seclark/FITSHandling/commit/f04a6e54c6624741e4f3077ba8ba96af620871ac

    for lb masks
    '''
    obj = coord(ras, decs, unit="deg", frame="icrs")
    obj = obj.galactic

    ls = obj.l.degree
    bs = obj.b.degree

    return ls, bs


def circ_kern(diameter):
    '''
    Performs a circle-cut of given diameter on inkernel.
    Outkernel is 0 anywhere outside the window.
    taken from https://github.com/seclark/RHT/blob/master/rht.py

    for umask step
    '''
    assert diameter % 2
    r = diameter // 2   # int(np.floor(diameter / 2))
    mnvals = np.indices((diameter, diameter)) - r
    rads = np.hypot(mnvals[0], mnvals[1])
    return np.less_equal(rads, r).astype(np.int)


def umask(data, radius=15, smr_mask=None):
    '''
    unsharp masking of a data slice
    taken from https://github.com/seclark/RHT/blob/master/rht.py
    ^conversion to binary data step skipped

    radius is set to 15 by default for GALFA data: diameter of 30 arcmin
    '''
    assert data.ndim == 2

    kernel = circ_kern(2 * radius + 1)
    outdata = scipy.ndimage.filters.correlate(data, kernel)

    # Correlation is the same as convolution here because kernel is symmetric
    # Our convolution has scaled outdata by sum(kernel), so we will divide out these weights.
    kernweight = np.sum(kernel)
    subtr_data = data - outdata / kernweight

    # No values < 0
    subtr_data[np.where(subtr_data < 0.0)] = 0

    # set NaNs to 0??
    # subtr_data[np.where(subtr_data == np.NaN)] = 0

    if smr_mask is None:
        return subtr_data
    else:
        return np.logical_and(smr_mask, subtr_data)


def add_tree_to_dict(tree, dictionary):
    '''
    Saves a tree to the dictionary and avoid name overlaps by adding the
    initial velocity and appending a letter if necessary

    Arguments:
        tree {[type]} -- [description]
        dict {[type]} -- [description]
    '''
    key = str(tree.getTreeMaskedArea2D())
    key += '_' + str(tree.getTreeStartingVelocity()) + '_0'

    while key in dictionary:
        key = tree_key_hash(key)

    dictionary[key] = tree


def tree_key_hash(original_key):
    '''
    hash them keys

    Arguments:
        original_key {[type]} -- [description]
    '''
    key_base = original_key.rsplit('_', 1)[0]
    key_num = int(original_key.rsplit('_', 1)[1])

    key_num += 1

    new_key = key_base + key_num
    return new_key


def delete_small_dead_trees(all_trees, length_cutoff=1, verbose=False):
    '''
    delete all the trees that have length = 1 and have ended

    Arguments:
        tree_dict {[type]} -- [description]
    '''
    bad_trees_keys = []

    for k in all_trees:
        if all_trees[k].has_ended and all_trees[k].length == length_cutoff:
            bad_trees_keys.append(k)
            continue

    if verbose:
        print "\t\t deleting %d trees that are small & dead" % len(bad_trees_keys)

    for k in bad_trees_keys:
        del all_trees[k]
