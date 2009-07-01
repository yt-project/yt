subroutine create_tree()
    use kdtree2_module
    use fKD_module
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    ! create a kd tree object

     tree2 => kdtree2_create(pos,sort=sort,rearrange=rearrange)  ! this is how you create a tree. 
     return

end subroutine create_tree


subroutine find_nn_nearest_neighbors()
     use kdtree2_module
     use fKD_module
     use kdtree2module
     use tree_nodemodule
     use intervalmodule

     integer :: k
     type(kdtree2_result),allocatable :: results(:) ! nearest neighbors
     !integer, parameter  :: nn ! number of nearest neighbors found


     allocate(results(nn)) 

     call kdtree2_n_nearest(tp=tree2,qv=qv,nn=nn,results=results) 

     dist = results%dis
     tags = results%idx

  !do k=1,nn
  !   print *, "k = ", k, " idx = ", tags(k)," dis = ", dist(k)
  !   print *, "x y z", pos(1,results(k)%idx), pos(2,results(k)%idx), pos(3,results(k)%idx)
  !enddo


     deallocate(results)
     return

end subroutine find_nn_nearest_neighbors

subroutine free_tree()
    use kdtree2_module
    use fKD_module
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    ! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
    ! TO MAKE THE TREE.  
    call kdtree2_destroy(tree2)
    
    ! Free up the particle memory
    ! not sure this works, might need to del(fKD.pos) in python space?
    deallocate(pos)
end subroutine free_tree

