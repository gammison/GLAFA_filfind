'''
This takes in a GALFA cube and processes it into a dict of dict
'''

# imports
import sys
import pickle
from astropy.io import fits
import mask_obj_node as maskNode
fil_finder_dir = '/Users/larryli/Documents/CC/16-17/research/GALFA_filfind/fil_finder'
sys.path.append(fil_finder_dir)
import filfind_class as filfind


def process_cube_filfind_struct(file_dir, file_name, v_range, x_range, y_range,
                                save_struct=True, verbose_process=False):
    '''
    Takes a GALFA data cube file, cuts it to the specified dimentions, and
    processes it slice by slice to find strucutres on each v slice with filfind.
    Each mask found by filfind on a single v slice is put into a dict with its
    masked_area_size as key and each dict is then put into a dict with its v
    index as key. The overall structure is either returned or pickled.
    '''

    # import cube
    cube_dir = file_dir
    cube_name = file_name
    full_cube, hdr = fits.getdata(cube_dir + cube_name, header=True)

    # full cube dimentions
    full_cube_shape = full_cube.shape
    full_v_channel_count = full_cube_shape[0]
    full_y_pixel_count = full_cube_shape[1]
    full_x_pixel_count = full_cube_shape[2]

    print "\n\tThere are %d velocity channels in total" % full_v_channel_count
    print "\n\tThe full image is %d by %d pixels" % (full_x_pixel_count, full_y_pixel_count)

    # cut cube based on provided x&y dimentions
    cut_cube = full_cube[:, y_range[0]:y_range[1], x_range[0]:x_range[1]]

    # run though the slices in v_range and find masks
    # store masks in nodes, and all nodes in a v slice in an dict by their masked_area_size
    # store that list in dict with v as key
    nodes_by_v_slice = {}
    for v in xrange(v_range[0], v_range[1]):
        v_slice = cut_cube[v, :, :]
        nodes_in_slice = {}
        print "\n\n\tworking on velocity slice %d" % v

        # puts slice into filfind
        fils = filfind.fil_finder_2D(v_slice, header=hdr, beamwidth=10.0, glob_thresh=20,
                                     distance=100, flatten_thresh=95, standard_width=1.1,
                                     size_thresh=1000)
        # note size_thresh, adapt_thresh, smooth_size, fill_hole_size can all be set by args
        mask_objs = fils.create_mask(verbose=verbose_process, regrid=False, border_masking=True,
                                     save_png=True, run_name=str(v), output_mask_objs=True)

        # put returned masks in a dict of mask_obj_nodes
        if mask_objs[0] == 0:
            print "\n\tNO objects in slice %d" % v
        else:
            for i in range(0, len(mask_objs[0])):
                this_mask_node = maskNode.MaskObjNode(mask_objs[0][i], mask_objs[1][i], v)
                nodes_in_slice[this_mask_node.masked_area_size] = this_mask_node

        # put that dict of mask_obj_nodes into a dict with v as key
        nodes_by_v_slice[v] = nodes_in_slice

    # exit options
    if save_struct:
        save_dir = '../pickled_dicts'
        save_name = cube_name.split('.')[0] + '_' + str(v_range) + '_' + str(x_range) + '_' + str(y_range)
        pickle.dump(nodes_by_v_slice, open(save_dir + save_name, 'wb'))

    else:
        return nodes_by_v_slice
