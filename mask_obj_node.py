import numpy as np


class MaskObjNode:
    '''
    This node object is made to contain the masks produced by fil_finder.
    Each node takes in a mask(square), the corners of the mask, and the v
    index, other fields are then calculted.

    It is important to note the differnce between the area of a mask and the
    masked area: since the masks are aways squares, the area of a mask is just
    calculated from dimentions of the corners. The masked area, however, is
    calculated by counting the amount of pixels that's masked (marked TRUE in
    the np array).

    '''

    def __init__(self, mask_obj, corners, v_slice_index):
        # mask is a 2d np bool array
        self.mask = mask_obj

        # corners are [(bottomleft)[y,x],(topright)[y,x]]
        # organize corners in to list of lists instead of list of tups
        # and convert into [x,y]
        self.corner_BL = [corners[0][1], corners[0][0]]
        self.corner_TR = [corners[1][1], corners[1][0]]
        self.corners = [self.corner_BL, self.corner_TR]

        self.v_slice_index = [v_slice_index, ]

        self.visited = False
        self.mask_size = self.checkAreaSize()
        self.masked_area_size = self.checkMaskedAreaSize()

    def mergeNode(self, other_node):
        '''
        Merge other_node with self. Append to self.v_slice_index if needed and
        create the combined OR mask of self.mask and other_node.mask. The
        combined or mask is set as the new self.mask and other attributes are
        fixed
        '''
        if self.v_slice_index[-1] != other_node.v_slice_index[0]:
            self.v_slice_index.append(other_node.v_slice_index)

        combined_or_mask = self.combineMask(other_node, merge_type='OR')
        combined_masked_area_size = self.checkMaskedAreaSize(combined_or_mask)
        new_corners = self.matchCorners(other_node)

        self.mask = combined_or_mask
        self.corners = new_corners
        self.corner_BL = new_corners[0]
        self.corner_TR = new_corners[1]
        self.mask_size = self.checkAreaSize(new_corners)
        self.masked_area_size = combined_masked_area_size

        return True

    def mergeNodeAlt(self, other_node):
        '''
        alt mergeNode function
        '''
        '''
        comboMask = [[0 for x in range(maxBR_X-TL_X)] for y in range(BR_Y - TL_Y)]
        radj_1 = TL_Y - mask1.leftcorner.y
        cadj_1 = TL_X - mask1.leftcorner.x

        for r in range(len(mask1)):
            for c in range(len(mask1[r])):
                comboMask[r-radj1][c-cadj1] = mask1[r][c]

        radj_2 = TL_Y - mask2.leftcorner.y
        cadj_2 = TL_X - mask2.leftcorner.x

        for r in range(len(mask2)):
            for c in range(len(mask2[r])):
                if mask2[r][c] == 1:
                    comboMask[r-radj2][c-cadj2] = mask2[r][c]

        return new mask(comboMask,TL_X,TL_Y,BL_X,BL_Y)
        '''

    def checkMaskOverlap(self, other_node, overlap_thresh):
        '''
        First check if the coners of self.mask and other_node.mask overlap, if
        the coners do overlap (meaning some part of the two squares overlap) we
        then check if actual masks overlap. To check for actual overlaps a
        combined AND mask is first made, and the masked area calculated. The
        masked area of the combined AND mask is then compared to the masked
        area of the imput masks. If the overlap is greater than overlap_thresh
        (# of pixels) then return True
        '''
        if self.checkCornersOverlap(other_node) == False:
            return False

        combined_and_mask = self.combineMask(other_node, merge_type='AND')
        combined_masked_area_size = float(self.checkMaskedAreaSize(combined_and_mask))

        if combined_masked_area_size / self.masked_area_size >= overlap_thresh:
            return True
        elif combined_masked_area_size / other_node.masked_area_size >= overlap_thresh:
            return True
        else:
            return False

    def combineMask(self, other_node, merge_type='AND'):
        '''
        Merge self.mask and other_node.mask. The corners of the 2 masks are
        first matched to fined the smallest square that contains both masks,
        a new mask is then made to that dimention. Elements of the new mask are
        filled basked on the merge_type ('AND' or 'OR').

        Return the new combined mask
        '''
        new_corners = self.matchCorners(other_node)
        r_dim = new_corners[1][1] - new_corners[0][1]
        c_dim = new_corners[1][0] - new_corners[0][0]

        combined_mask = np.zeros((r_dim, c_dim), dtype=bool)

        expanded_self_mask = self.expandMask(new_corners)
        expanded_other_mask = other_node.expandMask(new_corners)

        if merge_type == 'AND':
            combined_mask = np.bitwise_and(expanded_self_mask, expanded_other_mask)
        elif merge_type == 'OR':
            combined_mask = np.bitwise_or(expanded_self_mask, expanded_other_mask)

        return combined_mask

    def checkCornersOverlap(self, other_node):
        '''
        Check if self.mask and other_node.mask have any overlap
        '''
        if other_node.corner_BL[0] >= self.corner_TR[0] or \
           other_node.corner_TR[0] <= self.corner_BL[0]:
            return False
        if other_node.corner_BL[1] >= self.corner_TR[1] or \
           other_node.corner_TR[1] <= self.corner_BL[1]:
            return False

        return True

    def matchCorners(self, other_node):
        '''
        Compare the self.mask and other_node.mask and pick out the smallest
        square that contain both masks.

        Return the corners of that square.
        '''
        BL_X = min(self.corner_BL[0], other_node.corner_BL[0])
        BL_Y = min(self.corner_BL[1], other_node.corner_BL[1])

        TR_X = max(self.corner_TR[0], other_node.corner_TR[0])
        TR_Y = max(self.corner_TR[1], other_node.corner_TR[1])

        return [[BL_X, BL_Y], [TR_X, TR_Y]]

    def checkMaskedAreaSize(self, mask=None):
        '''
        Calculate the amount of pixels that are masked by self.mask. If a new
        mask is provided then that mask is used

        Return the # of pixels masked (as opposed to the area of the mask)
        '''
        if mask == None:
            mask = self.mask

        return np.size(np.where(mask == True)[0])

    def expandMask(self, new_corners):
        '''
        Expand self.mask so its corners match the new_corners provided. Masks
        are paded with 0s with the numpy function pad()
        Return new expanded mask
        '''
        mask = self.mask
        old_corners = self.corners

        t_Pad = abs(new_corners[1][1] - old_corners[1][1])
        b_Pad = abs(new_corners[0][1] - old_corners[0][1])
        l_Pad = abs(new_corners[0][0] - old_corners[0][0])
        r_Pad = abs(new_corners[1][0] - old_corners[1][0])

        return np.lib.pad(mask, ((b_Pad, t_Pad), (l_Pad, r_Pad)), 'constant', constant_values=0)

    def checkAreaSize(self, corners=None):
        '''
        Calculate the area of self.mask by looking at its dimentions from its
        corners. If corners are provided then an area based on that corner is
        calculated.

        Return the area of the mask (in pixels^2)
        '''
        if corners == None:
            corners = self.corners

        return (corners[1][0] - corners[0][0]) * (corners[1][1] - corners[0][1])
