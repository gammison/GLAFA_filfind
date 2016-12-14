import numpy as np

class MaskObjNode:
	'''
	'''

	def __init__(self, mask_obj, corners, v_slice_index):
		# mask is a 2d np bool array
		self.mask = mask_obj
		
		# corners are [(bottomleft)[y,x],(topright)[y,x]]
		# organize corners in to list of lists instead of list of tups
		self.corner_BL = [corners[0][1], corners[0][0]]
		self.corner_TR = [corners[1][1], corners[1][0]]
		self.corners = [self.corner_BL, self.corner_TR]

		self.v_slice_index = [v_slice_index,]

		self.visited = False
		self.mask_size = self.checkAreaSize()
		self.masked_area_size = self.checkMaskedAreaSize()



	def mergeNode(self, other_node):
		if self.v_slice_index[-1] != other_node.v_slice_index[0]:
			self.v_slice_index.append(other_node.v_slice_index)

		combined_or_mask = self.combineMask(other_node, merge_type='OR')
		combined_masked_area_size = checkMaskedAreaSize(combined_or_mask)
		new_corners = self.matchCorners(other_node)

		self.mask = combined_or_mask
		self.corners = new_corners
		self.corner_BL = new_corners[0]
		self.corner_TR = new_corners[1]
		self.mask_size = checkAreaSize(new_corners)
		self.masked_area_size = combined_masked_area_size

		return True


	def mergeNodeAlt(self, other_node):
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
		if self.checkCornersOverlap(other_node) == False:
			return False

		combined_and_mask = self.combineMask(other_node, merge_type='AND')
		combined_masked_area_size = float(checkMaskedAreaSize(combined_and_mask))

		if combined_masked_area_size/self.masked_area_size >= overlap_thresh:
			return True
		elif combined_masked_area_size/other_node.masked_area_size >= overlap_thresh:
			return True
		else:
			return False


	def combineMask(self, other_node, merge_type='AND'):
		new_corners = self.matchCorners(other_node)
		r_dim = new_corners[1][1] - new_corners[0][1]
		c_dim = new_corners[1][0] - new_corners[0][0]

		combined_mask = np.zeros((r_dim, c_dim), dtype=bool)

		expanded_self_mask = self.expandMask(self.corners, new_corners)
		expanded_other_mask = other_node.expandMask(other_node.corners, new_corners)

		if merge_type == 'AND':
			combined_mask = np.bitwise_and(expanded_self_mask, expanded_other_mask)
		elif merge_type == 'OR':
			combined_mask = np.bitwise_or(expanded_self_mask, expanded_other_mask)

		return combined_mask


	def checkCornersOverlap(self, other_node):
		if other_node.corner_BL[0] >= self.corner_TR[0] or \
			other_node.corner_TR[0] <= self.corner_BL[0]:
			return False
		if other_node.corner_BL[1] >= self.corner_TR[1] or \
			other_node.corner_TR[1] <= self.corner_BL[1]:
			return False

		return True


	def matchCorners(self, other_node):
		BL_X = min(self.corner_BL[0], other_node.corner_BL[0])
		BL_Y = min(self.corner_BL[1], other_node.corner_BL[1])

		TR_X = min(self.corner_TR[0], other_node.corner_TR[0])
		TR_Y = min(self.corner_TR[1], other_node.corner_TR[1])

		return [[BL_X, BL_Y], [TR_X, TR_Y]]


	def checkMaskedAreaSize(self, mask=None):
		if mask == None:
			mask = self.mask

		return np.size(np.where(mask == True)[0])


	def checkAreaSize(self, corners=None):
		if corners == None:
			corners = self.corners

		return (corners[1][0]-corners[0][0]) * (corners[1][1]-corners[0][1])


	def expandMask(self, new_corners):
		mask = self.mask
		old_corners = self.corners

		t_Pad = abs(new_corners[1][1] - old_corners[1][1])
		b_Pad = abs(new_corners[0][1] - old_corners[0][1])
		l_Pad = abs(new_corners[0][0] - old_corners[0][0])
		r_Pad = abs(new_corners[1][0] - old_corners[1][0])

		return np.lib.pad(mask, ((b_Pad,t_Pad),(l_Pad,r_Pad)), 'constant', constant_values=0)

