import numpy as np
from unyt.array import ustack

from yt.units import YTArray
from yt.utilities.amr_kdtree.api import AMRKDTree


def _accumulate_vector_field(path, field_vals):
    r"""
    Integrates the given vector field along the given path p.

    The integral is done in a piecewise manner (segment-by-segment) so as
    to be able to store the accumulated values from each of the previous
    segments.

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

    Since the integral involves evaluating the intervals :math:`\Delta x, \Delta y, \Delta z`, we have to decide which endpoint of each segment contributes its field
    value to the calculation. Here we use the average of the two values.

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
    f = (field_vals[1:] + field_vals[:-1]) / 2.0
    accum = np.cumsum(np.diag(np.dot(f.d, (path[1:].d - path[:-1].d).T)))
    accum = YTArray(accum, field_vals.units * path.units)
    return accum


def _accumulate_scalar_field(p, field_vals):
    r"""
    This function integrates a scalar field along a path:

    ..math::

        I = \int_C \phi(x1,x2,...,xn)d\vec{r}

    With this, the result will be a vector. As such, we can take advantage
    of numpy broadcasting to vectorize this calculation. If `p` is (N,d),
    where N is the number of points parameterizing the path and d is the
    number of spatial dimensions and phi is (N, 1) then the result will
    be (N,d) and have the form

    ..math::

        \begin{pmatrix}
        \phi_1 r_{1,x} & \phi_1 r_{1,y} & \phi_1 r_{1,z} \\
        \phi_2 r_{2,x} & \phi_2 r_{2,y} & \phi_2 r_{2,z} \\
        \phi_3 r_{3,x} & \phi_3 r_{3,y} & \phi_3 r_{3,z} \\
        \vdots & \ddots & \vdots
        \end{pmatrix}

    Where each column corresponds to the components of the resulting vector.
    We then do the cumulative sum along these columns to get the accumulation
    of the integral.

    As with `_accumulate_vector_field`, since we are dealing with intervals,
    we use the average value of the field determined from the values at the
    endpoints of each segment.

    Parameters
    ----------
    p : YTArray
        The path to be integrated along

    field_vals : YTArray
        An array containing the values of the scalar field to be
        integrated at the location of the starting point of each
        path segment as well as the endpoint for the last segment

    Returns
    -------
    accum : YTArray
        The cumulative value of the field integral at each path
        segment
    """
    f = (field_vals[1:] + field_vals[:-1]) / 2.0
    accum = np.cumsum(f.d * (p[1:].d - p[:-1].d), axis=0)
    return YTArray(accum, field_vals.units * p.units)


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
    """

    def __init__(self, paths, ds):
        self.ds = ds
        self.paths = self._verify_paths(paths)
        self.ad = self.ds.all_data()
        self.left_edge = self.ds.domain_left_edge
        self.right_edge = self.ds.domain_right_edge
        self.accum = []

    def _verify_paths(self, paths):
        """
        The ``Accumulators`` object can accept either a list of paths or
        a single path to integrate over. Here we make sure that
        ``paths`` is one of those two options and that the given paths
        are of the correct shape (either two or three dimensions,
        depending on the dimensionality of the given dataset).
        """
        verified_paths = []
        # Case 1: paths is a list
        if isinstance(paths, list):
            try:
                assert len(paths) > 0
            except AssertionError:
                raise ValueError("List of paths is empty.")
            for p in paths:
                verified_paths.append(self._verify_path(p))
        # Case 2: paths is an array
        elif isinstance(paths, np.ndarray):
            verified_paths.append(self._verify_path(paths))
        else:
            raise TypeError("Paths must be either an array or list of arrays.")
        return verified_paths

    def _verify_path(self, p):
        if isinstance(p, np.ndarray):
            # Make sure the shape is 2D
            try:
                assert len(p.shape) == 2
            except AssertionError:
                raise ValueError("Path shape should be (N, ds.dimensionality)")
            # Make sure the dimensionality of path points matches ds
            try:
                assert p.shape[1] == self.ds.dimensionality
            except AssertionError:
                raise ValueError("Point and ds dimensions don't match.")
            # Make sure there are at least two points on the path
            try:
                assert p.shape[0] > 1
            except AssertionError:
                raise ValueError("Need at least two points on path.")
            # Check units
            if hasattr(p, "units") and not p.units.is_dimensionless:
                return p.to(self.ds.length_unit)
            else:
                return self.ds.arr(p, "code_length")
        else:
            raise TypeError("Path must be either a ndarray or YTArray")

    def _get_tree(self, field):
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
        dims = data[0].shape
        # Put the data in a more convenient form: dims x ndims array
        data = [d.flatten() for d in data]
        data = ustack(data, axis=1)
        # Cell width
        node_left_edge = self.ds.arr(node.get_left_edge(), "code_length")
        node_right_edge = self.ds.arr(node.get_right_edge(), "code_length")
        cell_size = (node_right_edge - node_left_edge) / dims
        while idx < npts:
            # Make sure point is within domain
            left_check = path[idx] < self.left_edge
            right_check = path[idx] >= self.right_edge
            if np.any(left_check | right_check):
                msg = f"Point `{path[idx]}` at index `{idx}` outside domain bounds."
                msg += f"LE: `{self.left_edge}`, RE: `{self.right_edge}`"
                raise ValueError(msg)
            # Figure out which cell in the node the point falls within
            # Origin of node can be offset from origin of origin of volume,
            # so we have to subtract it off to get the right cell indices
            cell_ind = ((path[idx] - node_left_edge) / cell_size).astype("i8")
            # Access the value of the field at that index. Accessing a single
            # element of a multi-dimensional array using another array is
            # problematic. Flatten and use a row-major index, indstead
            I = np.ravel_multi_index(cell_ind, dims)
            # Get the value of each component for the current cell
            vals[idx] = data[I]
            # See if next point is still within the same node (if the next point
            # is still in range of the array)
            idx += 1
            if idx != npts:
                left_check = path[idx] < node_left_edge
                right_check = path[idx] >= node_right_edge
                if np.any(left_check | right_check):
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
        if is_vector is None:
            raise ValueError("is_vector parameter not set.")
        if not isinstance(field, list):
            field = [
                field,
            ]
        self.accum = []
        # Build tree and add field to it
        tree = self._get_tree(field)
        # Loop over each path the field is to be accumulated along
        for p in self.paths:
            # Get the value of the field at each point on the path
            npts = len(p)
            if npts < 2:
                raise ValueError("Invalid path. Must have at least 2 points.")
            if is_vector:
                values = YTArray(np.zeros(p.shape), self.ad[field[0]].units)
            else:
                values = YTArray(np.zeros((npts, 1)), self.ad[field[0]].units)
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
