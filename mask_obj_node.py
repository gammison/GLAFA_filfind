import numpy as np

class MaskObjNode:

	def __init__(self, mask_obj, corners, v_slice_index):
		self.mask_obj = mask_obj #2d array
		self.corners = corners
		self.v_slice_index = v_slice_index

		self.visited = False
		self.mask_size = (np.where(mask_obj == True)).size


	def mergeNode(self, other):
		pass

	def hasOverlap():
		pass