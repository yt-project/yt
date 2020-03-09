import numpy as np
import unyt
from unyt.array import ustack

from yt.funcs import ensure_list
from yt.units.yt_array import YTArray
from yt.utilities.amr_kdtree.api import AMRKDTree


def _accumulate_vector_field(path, field_vals):
    r"""
    This function integrates the given vector field along the given
    path p. The integral is done in a piecewise manner
    (segment-by-segment) so as to be able to store the accumulated
    values from each of the previous segments.

    For a vector field, the line integral is:

    ..math::

        I = \int_C \vec{a} \cdot d\vec{r}

    where :math:`C` is the path being integrated along, :math:`\vec{a}` is
    the vector field being integrated, and :math:`d\vec{r}` points along the
    path :math:`C`. This is equivalent to:
    
    ..math::

        \lim_{N \rightarrow \infty} \sum_{p=1}^N \vec{a}(x_p, y_p, ...) \cdot \Delta \vec{r}_p
        
    with the understanding that :math:`\Delta \vec{r}_p \rightarrow 0` as
    :math:`N \rightarrow \infty`.
        
    The vector pointing along the segment connecting two adjacent points is:
    
    ..math::
    
        \Delta\vec{r}_p = (\vec{r}_{p+1} - \vec{r}_p)
        
    so the full dot product can be written out as:

    ..math::
        I = \lim_{N \rightarrow \infty}\sum_{p=1}^N\sum_{i=1}^n a_{p,i}(r_{p+1,i}-r_{p,i})
        
    where :math:`n` is the number of dimensions.

    This can be done as a matrix operation. If :math:`\vec{r} \equiv \vec{r}_{p+1} - \vec{r}_p`
    then the above sum is (leaving out the limit):
    
    ..math::
    
        I = \vec{a}_1 \dot \vec{r}_1 + \vec{a}_2 \cdot \vec{r}_2 \ldots \vec{a}_N \cdot \vec{r}_N
        
    If we write the matrix :math:`A = ([a1x, a1y, ...], [a2x, a2y, ...], ...)`
    and :math:`R = ([r1x, r1y, ...], [r2x, r2y, ...], ...)^T`, then the dot
    products in the sum are the diagonal elements of the resulting matrix
    multiplication :math:`AR`. The accumulation is then obtained by doing a
    cumsum of this diagonal.

    Parameters
    ----------
    p : YTArray 
        The path to be integrated along

    field_vals : YTArray 
        An array containing the components of the vector field to be
        integrated. The values are sampled at the starting point of
        each path segment as well as the endpoint for the last segment

    Returns
    -------
        accum : YTArray 
            The cumulative value of the field integral at each path
            segment
    """
    # https://en.wikipedia.org/wiki/Line_integral
    # np.dot doesn't combine units correctly and there does not appear to be a
    # unyt implementation of dot, so units will be handled manually
    # for the time being
    accum = np.cumsum(np.diag(np.dot(field_vals[:-1].d, (path[1:].d - path[:-1].d).T)))
    accum = YTArray(accum, field_vals.units * path.units)
    return accum


def _accumulate_scalar_field(p, field_vals):
    r"""
    This function integrates a scalar field along a path. It uses a
    similar method to that in _accumulate_vector_field, but the integral
    is now:
    
    ..math::
    
        I = \int_C \phi(x1,x2,...,xn)d\vec{r}

    Parameters
    ----------
    p : YTArray 
        The path to be integrated along

    fieldVals : YTArray 
        An array containing the values of the scalar field to be
        integrated at the location of the starting point of each
        path segment as well as the endpoint for the last segment

    Returns
    -------
    accum : YTArray 
        The cumulative value of the field integral at each path
        segment 
    """
    # https://en.wikipedia.org/wiki/Line_integral
    # np.linalg.norm returns a ndarray, so when multiplied by field_vals, which
    # is a YTArray, this leads to incorrect units. There does not appear to be
    # a unyt implementation of norm, as far as I'm aware, so units will be
    # handled manually for the time being
    accum = np.cumsum(field_vals[:-1].d * np.linalg.norm(p[1:].d - p[:-1].d, axis=1))
    accum = YTArray(accum, field_vals.units * p.units)
    return accum


def get_row_major_index(ncells, cell_ind):
    """
    Converts the cell indices (i_1, i_2, ..., i_N) to the corresponding
    row-major index I.

    The field data for each cell in the node are stored in an N dimensional
    array with shape (n_1, n_2, ..., n_N), where n_i is the number of cells along
    dimension i. The N dimensional index of the current cell under consideration,
    `cell_ind` is (i_1, i_2, ..., i_N). However, accessing the value of a single
    element in a multi-dimensional array (unflattened data) with another array
    (cell_ind) is problematic.

    As such, we flatten the data, which puts the data into row-major format.
    The current cell under consideration can then be accessed by converting
    cell_ind to the corresponding row-major index I.

    The location I in the row-major array of the cell given by (1_1, i_2, ..., i_N)
    is:

        I = i_N + n_N * (i_{N-1} + n_{N-1} * (...(i_2 + n_2 * i)...)

    which in three dimensions is:

        I = k + n_z * (j + n_y * i)

    Parameters
    ----------
    ncells : tuple
        Contains the number of cells along each dimension.

    cell_ind : tuple
        Contains the integer values of the index of the current cell in each
        dimension (i.e., this is (i_1, i_2, ..., i_N)).

    Returns
    -------
    I : int
        The row-major index corresponding to cell_ind.
    """
    I = cell_ind[0]
    for idx in range(1, len(ncells)):
        I = cell_ind[idx] + ncells[idx] * I
    return I


class Accumulators:
    r"""
    Container for creating and storing the path integrals of various
    fields along various paths.

    The class takes in a list of user-defined paths and a dataset. From
    these, the user can compute the integral of any field in the dataset
    along the given paths. The results of these path integrals are
    stored in a cumulative fashion, which means that, once computed, the
    accumulated value of the field is stored at each point along the
    path, allowing the user to query this information.

    Attributes
    ----------
    pass

    Methods
    -------
    pass

    Examples
    --------
    pass
    """
    def __init__(self, paths, ds):
        self.paths      = []
        self.ds         = ds
        self.ad         = ds.all_data()
        self.accum      = []
        self.left_edge  = self.ds.domain_left_edge
        self.right_edge = self.ds.domain_right_edge
        # Make sure that the path has proper units. If no units are
        # specified, assume code units
        for p in paths:
            if not isinstance(p, YTArray) or p.units == unyt.dimensionless:
                self.paths.append(self.ds.arr(p, 'code_length'))
            else:
                self.paths.append(p.to(self.ds.length_unit))
                    
                
                

    def _get_tree(self, field, is_vector):
        r"""
        Creates an AMRKDTree and then adds the field to it.

        Parameters
        ----------
        field : list 
            Either the scalar field to add to the tree or an iterable
            containing the components of the vector field to add.

        is_vector : bool
            If True, then `field` is a vector and it is assumed that it's
            an iterable containing the x, y, and z components of the field.
            Otherwise, the field is assumed to be a scalar field.

        Returns
        -------
        tree : AMRKDTree
            An AMRKDTree for the dataset containing the given field values.
        """
        tree = AMRKDTree(self.ds)
        use_log = [False for i in range(len(field))]
        field = [self.ad._determine_fields(f)[0] for f in field]
        tree.set_fields(field, use_log, False)
        return tree

    def _get_path_field_values(self, tree, path, idx, vals, npts):
        r"""
        Determines the value of the field at each point along the path.

        Parameters
        ----------
        tree : AMRKDTree
            An AMRKDTree for the dataset containing the given field values.

        path : np.ndarray
            The N x 3 array containing the x, y, and z coordinates of the
            N points used to define the path.

        Returns
        -------
        path_field_values : np.ndarray
            The array containing the field values at each point along the
            path.
        """
        node = tree.locate_node(path[idx])
        # Data is a list, with one element for each field component
        data = node.data.my_data
        # Number of cells in the node
        ncells = data[0].shape
        # Put the data in a more convenient form: ncells x ndims array
        data = [d.flatten() for d in data]
        data = ustack(data, axis=1)
        # Cell width
        node_left_edge = self.ds.arr(node.get_left_edge(), 'code_length')
        node_right_edge = self.ds.arr(node.get_right_edge(), 'code_length')
        cell_size = (node_right_edge - node_left_edge) / ncells
        while idx < npts:
            # Make sure point is within domain
            left_check = path[idx] < self.left_edge
            right_check = path[idx] >= self.right_edge
            if np.sum(np.logical_or(left_check, right_check)):
                msg = f"Point `{path[idx]}` at index `{idx}` outside domain bounds."
                msg += f"LE: `{self.left_edge}`, RE: `{self.right_edge}`"
                raise ValueError(msg)
            # Figure out which cell in the node the point falls within
            # Origin of node can be offset from origin of origin of volume,
            # so we have to subtract it off to get the right cell indices
            cell_ind = ((path[idx] - node_left_edge) / cell_size).astype(int)
            # Access the value of the field at that index. Accessing a single
            # element of a multi-dimensional array using another array is
            # problematic. Flatten and use a row-major index, indstead  
            I = get_row_major_index(ncells, cell_ind)
            # Get the value of each component for the current cell
            vals[idx] = data[I] 
            # See if next point is still within the same node (if the next point
            # is still in range of the array)
            idx += 1
            if idx != npts:
                left_check = path[idx] < node_left_edge
                right_check = path[idx] >= node_right_edge
                if np.sum(np.logical_or(left_check, right_check)):
                    vals, idx = self._get_path_field_values(tree, path, idx, vals, npts)
        return vals, idx
                
    def accumulate(self, field, is_vector=None):
        r"""
        This function is the driver function for integrating the desired
        field along each path in the bundle.

        Parameters
        ----------
        field : field name, iterable of field names
            The field to integrate along the paths. If this is a vector
            field, field should be an iterable containing the components
            of the field (i.e, field[0] is the x-component data, field[1]
            is the y-component, etc.). This is the name(s) of the field
            (components), e.g., ('gas', 'density') for a scalar or
            [('gas', 'x'), ('gas', 'y'), ('gas', 'z')] for a vector.

        is_vector : bool
            If True, field is a vector field. If False, then field is
            a scalar field.

        Raises
        ------
        ValueError
            If is_vector is not set or if one of the givne pathes is
            too short (less than two points).
        """
        import pdb; pdb.set_trace()
        if is_vector is None:
            raise ValueError("`is_vector` parameter not set.")
        field = ensure_list(field)
        self.accum = []
        # Build tree and add field to it
        tree = self._get_tree(field, is_vector)
        # Loop over each path the field is to be accumulated along
        for p in self.paths:
            # Get the value of the field at each point on the path
            npts = len(p)
            if npts < 2:
                raise ValueError("Invalid path. Must have at least 2 points.")
            values = YTArray(np.zeros(p.shape), self.ad[field[0]].units) 
            values, _ = self._get_path_field_values(tree, p, 0, values, npts)
            if is_vector:
                self.accum.append(_accumulate_vector_field(p, values))
            else:
                self.accum.append(_accumulate_scalar_field(p, values))

    def clear(self):
        r"""
        Clears out the accum attribute.
        """
        self.accum = []

    def add_path(self, path):
        r"""
        Allows the user to add a path to the accumulator after instantiation.
        """
        self.paths.append(path)
