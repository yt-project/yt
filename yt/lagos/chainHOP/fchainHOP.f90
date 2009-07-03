! A fortran implementation of chainHOP
! Stephen Skory, LCA/CASS/UCSD 2009

! This is the main subroutine that calls other subroutines
subroutine fchainHOP()
    use kdtree2_module ! from the kdtree2 source file
    use fchainHOPmodule ! from the .v file
    use NN_module ! in this file
    use sub_module ! in this file
    
    ! ----- INPUTS ------
    ! (set in Python, kept in the fchainHOPmodule from the .v file)
    ! pos - real array(3,nparts), the positions of the particles
    ! num_neighbors - integer, number of neighbors to find in the nearest
    !   neighbor search, which also determines the number of particles to
    !   calculate the smoothing kernel over
    ! min_bounds, max_bounds - real array(3), the corners of the 'real' volume
    !   inspected by this instance of HOP. The particles outside of this box
    !   are in the padding.
    ! mass - real array(nparts), masses of the particles, normalized to 1
    ! threshold - real, threshold for group creation (delta_outer).
    ! nparts - integer, number of particles
    ! npadded - integer, maximum number of padded particles that can go into
    !   padded_particles
    ! rearrange - logical, whether the kdtree should rearrange the position
    !   array internally for greater speed at the cost of memory usage.
    ! sort - logical, whether the results of a nearest neighbor search should
    !   be returned sorted in ascending order of distance. On is slower, of
    !   course.
    ! (?) period - real array(3), the period of the periodicity of the volume
    !
    ! ------ OUTPUTS -----
    ! (instantiated in Python)
    ! tags - integer array(nparts), the chainIDs for each particle, including 
    !   boundary particles.
    ! dens - real array(nparts), the densities for each particle
    ! padded_particles - integer array(npadded), the particleID of particles
    !   in the padding that are part of chains.
    !
    
    ! ---------------- Let's go! --------
    
    integer :: i, a_count
    
    ! kd tree object
    type(kdtree2), pointer :: tree2
    ! query particle point
    real(kdkind), allocatable :: query_vec(:)
    ! nearest neighbors object
    type(kdtree2_result),allocatable :: results(:)
    type(typeNN), pointer :: NN

    allocate(results(num_neighbors)) 
    allocate(query_vec(3))    
    ! create the tree
    print *, "building the tree...."
    tree2 => kdtree2_create(pos,sort=sort,rearrange=rearrange)
   
    ! fill the NN object
    NN => init_NN()
    
    ! loop over the particles finding nearest neighbors, and the density
    ! at that point.
    print *, "finding nearest neighbors/density..."
    do i=1,nparts
        query_vec(:) = pos(:,i)
        call kdtree2_n_nearest(tp=tree2,qv=query_vec,nn=num_neighbors,&
            results=results)
        ! copy over the results
        NN%NNlist_ID(i,:) = results%idx
        NN%NNlist_d(i,:) = results%dis
        ! now find the density at this particle
        call smDensitySym(i,NN)
    end do
    ! done searching, clean up
    deallocate(results)
    deallocate(query_vec)
    call kdtree2_destroy(tree2)
    
    a_count = count(MASK = NN%density .GE. threshold)
    print *, "There are ", a_count, " particles above the threshold"
    
    ! Each particle has a density, a list of NN, let's find their densest
    ! nearest neighbor
    print *, "Finding densest nearest neighbors..."
    do i=1,nparts
        call densestNN(i,NN)
    end do
    
    a_count = 0
    do i=1,nparts
        if (NN%densestNN(i) .EQ. i) a_count = a_count + 1
    end do
    print *, "There are ", a_count, " self-densest particles."
    
    ! First round, build chains of links
    !print *, "Building initial chains..."
    !call build_chains()
    
end subroutine fchainHOP

