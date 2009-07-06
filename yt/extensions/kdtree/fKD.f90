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

subroutine find_all_nn_nearest_neighbors()
    ! for all particles in pos, find their nearest neighbors and return the
    ! indexes and distances as big arrays
    use kdtree2_module
    use fKD_module
    use kdtree2module
    use tree_nodemodule
    use intervalmodule

    integer :: k
    type(kdtree2_result),allocatable :: results(:) ! nearest neighbors
    allocate(results(nn))
    
    do k=1,nparts
        qv(:) = pos(:,k)
        call kdtree2_n_nearest(tp=tree2,qv=qv,nn=nn,results=results)
        nn_dist(:,k) = results%dis
        nn_tags(:,k) = results%idx
    end do
    
    deallocate(results)
    return

end subroutine find_all_nn_nearest_neighbors

subroutine chainHOP_tags_dens()
    ! for all particles in pos, find their nearest neighbors, and calculate
    ! their density.
    use kdtree2_module
    use fKD_module
    use kdtree2module
    use tree_nodemodule
    use intervalmodule

    integer :: k, pj, i
    real :: ih2, fNorm, r2, rs
    type(kdtree2_result),allocatable :: results(:) ! nearest neighbors
    allocate(results(nn))
    
    do k=1,nparts
        qv(:) = pos(:,k)
        call kdtree2_n_nearest(tp=tree2,qv=qv,nn=nn,results=results)
        nn_tags(:,k) = results%idx
        nn_dist(:,k) = results%dis
        
        ! calculate the density for this particle
        ih2 = 4.0/nn_dist(nn,k)
        fNorm = 0.5*sqrt(ih2)*ih2/3.1415926535897931
        do i=1,nn
            pj = nn_tags(i,k)
            r2 = nn_dist(i,k)*ih2
            rs = 2.0 - sqrt(r2)
            if (r2 < 1.0) then
                rs = (1.0 - 0.75*rs*r2)
            else
                rs = 0.25*rs*rs*rs
            end if
            rs = rs * fNorm
            dens(k) = dens(k) + rs * mass(pj)
            dens(pj) = dens(pj) + rs * mass(k)
        end do
    
    end do
    
    deallocate(results)
    return

end subroutine chainHOP_tags_dens

subroutine free_tree()
    use kdtree2_module
    use fKD_module
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    ! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
    ! TO MAKE THE TREE.  
    call kdtree2_destroy(tree2)
    
    ! The data to make the tree has to be deleted in python BEFORE calling
    ! this!
end subroutine free_tree

