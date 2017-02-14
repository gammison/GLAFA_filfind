import numpy as np
import mask_obj_node

test1 = np.array([[False,False,False,False,False,False],
		 [True,True,True,False,False,False],
		 [True,True,True,False,False,False],
		 [False,False,False,False,False,False]])

test2 = np.array([[True,True,True],
		  [False,False,False]])

test3 = np.array([[False,False,False,False,False,False],
		 [True,True,True,False,True,False],
		 [True,True,True,True,True,False],
		 [False,False,False,False,False,False]])

test1_corner = [[0,0],[3,5]]
test2_corner = [[1,0],[2,2]]
test3_corner = [[0,0],[3,5]]

test1_vel = 1
test2_vel = 2
test3_vel = 3

node_1 = mask_obj_node.MaskObjNode(test1,test1_corner,test1_vel)
node_2 = mask_obj_node.MaskObjNode(test2,test2_corner,test2_vel)
node_3 = mask_obj_node.MaskObjNode(test3,test3_corner,test3_vel)

#print(node_2.expandMask([[0,0],[7,9]]))
#print(node_1.expandMask([[0,0],[7,9]]))
#print(node_2.matchCorners(node_1))
#combined_mask = node_1.combineMask(node_2,merge_type='AND')

#print(combined_mask)
print(node_1.masked_area_size)
print(node_2.masked_area_size)
#print(node_1.checkMaskOverlap(node_2,0.75))



