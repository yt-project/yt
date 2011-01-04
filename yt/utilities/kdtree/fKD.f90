subroutine create_tree(treeID)
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    integer, intent (in), optional         :: treeID
    
    integer :: ID
    
    if (present(treeID)) then
       ID = treeID
    else
       ID = -1
    end if
    
    ! create a kd tree object
    
    if (ID == 1) then
        t1%tree2 => kdtree2_create(t1%pos,sort=t1%sort,rearrange=t1%rearrange)
    elseif (ID == 2) then
        t2%tree2 => kdtree2_create(t2%pos,sort=t2%sort,rearrange=t2%rearrange)
    else
        tree2 => kdtree2_create(pos,sort=sort,rearrange=rearrange)
    end if

    return

end subroutine create_tree

subroutine add_tree(treeID)
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    integer :: treeID
    
    if (treeID == 1) then
        t1 => Newtree_set()
    elseif (treeID == 2) then
        t2 => Newtree_set()
    end if
    
    return

end subroutine add_tree

subroutine find_nn_nearest_neighbors()
     use kdtree2_module
     use fKD_module
     use tree_setmodule
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

subroutine find_many_nn_nearest_neighbors()
    ! Given an input array of query vectors (qv_many), find their
    ! nearest neighbors.
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    integer :: k, number
    type(kdtree2_result),allocatable :: results(:)

    allocate(results(nn))

    number = size(qv_many,2)

    do k=1, number
        qv(:) = qv_many(:,k)
        call kdtree2_n_nearest(tp=tree2,qv=qv,nn=nn,results=results)
        nn_tags(:, k) = results%idx
    end do
    
    deallocate(results)
    return

end subroutine find_many_nn_nearest_neighbors

subroutine find_r_nearest()
    ! Given an input array of a point (qv), find the number, fortran indices and
    ! squared distances to all neighbors within (radius) of (qv).
    ! The number of neighbors with indices and radius returned is set by
    ! (radius_n), and if this is smaller than (nfound) when returned,
    ! there are neighbors missing from the results.
    ! Store the results in (dist) and (tags).
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule
    
    integer :: k
    type(kdtree2_result),allocatable :: results(:) ! nearest neighbors
    
    allocate(results(radius_n))
    nfound = 0
    
    ! do it.
    call kdtree2_r_nearest(tp=tree2,qv=qv,r2=radius,nfound=nfound,nalloc=radius_n,results=results)
    
    tags(:) = results%idx
    dist(:) = results%dis
    
    deallocate(results)
    return

end subroutine find_r_nearest

subroutine find_many_r_nearest(treeID)
    ! Given an input array of a points (qv_many), find the number, fortran indices and
    ! squared distances to all neighbors within (radius) of (qv).
    ! The number of neighbors with indices and radius returned is set by
    ! (radius_n), and if this is smaller than (nfound) when returned,
    ! there are neighbors missing from the results.
    ! Store the results in (dist) and (tags).
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule

    integer, intent (in), optional         :: treeID
    
    integer :: ID
    integer :: k, number
    type(kdtree2_result),allocatable :: results(:) ! nearest neighbors

    
    if (present(treeID)) then
       ID = treeID
    else
       ID = -1
    end if    

    if (ID==1) then
        allocate(results(t1%radius_n))
        number = size(t1%qv_many,2)
        do k=1, number
            t1%qv(:) = t1%qv_many(:,k)
            t1%nfound = t1%nfound_many(k)
            call kdtree2_r_nearest(tp=t1%tree2,qv=t1%qv,r2=t1%radius,nfound=t1%nfound,nalloc=t1%radius_n,results=results)
            t1%nfound_many(k) = t1%nfound
            t1%nn_tags(:, k) = results%idx
            t1%nn_dist(:, k) = results%dis
        end do
    elseif (ID==2) then
        allocate(results(t2%radius_n))
        number = size(t2%qv_many,2)
        do k=1, number
            t2%qv(:) = t2%qv_many(:,k)
            t2%nfound = t2%nfound_many(k)
            call kdtree2_r_nearest(tp=t2%tree2,qv=t2%qv,r2=t2%radius,nfound=t2%nfound,nalloc=t2%radius_n,results=results)
            t2%nfound_many(k) = t2%nfound
            t2%nn_tags(:, k) = results%idx
            t2%nn_dist(:, k) = results%dis
        end do
    else
        allocate(results(radius_n))
        number = size(qv_many,2)
        do k=1, number
            qv(:) = qv_many(:,k)
            nfound = nfound_many(k)
            call kdtree2_r_nearest(tp=tree2,qv=qv,r2=radius,nfound=nfound,nalloc=radius_n,results=results)
            nfound_many(k) = nfound
            nn_tags(:, k) = results%idx
            nn_dist(:, k) = results%dis
        end do
    end if
    
    deallocate(results)
    return

end subroutine find_many_r_nearest

subroutine find_all_nn_nearest_neighbors()
    ! for all particles in pos, find their nearest neighbors and return the
    ! indexes and distances as big arrays
    use kdtree2_module
    use fKD_module
    use tree_setmodule
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

subroutine find_chunk_nearest_neighbors()
    ! for a chunk of the full number of particles, find their nearest neighbors
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule

    integer :: k
    type(kdtree2_result),allocatable :: results(:) ! nearest neighbors
    allocate(results(nn))
    do k=start,finish
        qv(:) = pos(:,k)
        call kdtree2_n_nearest(tp=tree2,qv=qv,nn=nn,results=results)
        chunk_tags(:,k - start + 1) = results%idx

    end do
    
    deallocate(results)
    return

end subroutine find_chunk_nearest_neighbors

subroutine chainHOP_tags_dens()
    ! for all particles in pos, find their nearest neighbors, and calculate
    ! their density. Return only nMerge nearest neighbors.
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule

    integer :: k, pj, i
    real :: ih2, fNorm, r2, rs
    integer, allocatable :: temp_tags(:)
    real, allocatable :: temp_dist(:)
    type(kdtree2_result),allocatable :: results(:) ! nearest neighbors
    allocate(results(nn))
    allocate(temp_tags(nn))
    allocate(temp_dist(nn))
    
    do k=1,nparts
        qv(:) = pos(:,k)
        
        call kdtree2_n_nearest(tp=tree2,qv=qv,nn=nn,results=results)
        temp_tags(:) = results%idx
        temp_dist(:) = results%dis
        
        ! calculate the density for this particle
        ih2 = 4.0/maxval(results%dis)
        fNorm = 0.5*sqrt(ih2)*ih2/3.1415926535897931
        do i=1,nn
            pj = temp_tags(i)
            r2 = temp_dist(i) * ih2
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

        ! record only nMerge nearest neighbors, but skip the first one which
        ! is always the self-same particle
        ! nn_tags(:,k) = temp_tags(2:nMerge)
    end do
    
    deallocate(results)
    deallocate(temp_dist)
    deallocate(temp_tags)
    return

end subroutine chainHOP_tags_dens

subroutine free_tree(treeID)
    use kdtree2_module
    use fKD_module
    use tree_setmodule
    use kdtree2module
    use tree_nodemodule
    use intervalmodule

    integer, intent (in), optional         :: treeID
    
    integer :: ID
    
    if (present(treeID)) then
       ID = treeID
    else
       ID = -1
    end if
    
    ! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
    ! TO MAKE THE TREE.  
    if (ID == 1) then
        call kdtree2_destroy(t1%tree2)
    elseif (ID == 2) then
        call kdtree2_destroy(t2%tree2)
    else
        call kdtree2_destroy(tree2)
    end if
    
    ! The data to make the tree has to be deleted in python BEFORE calling
    ! this!
end subroutine free_tree

