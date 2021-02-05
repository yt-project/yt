import numpy as np

from yt.funcs import mylog


def receive_and_reduce(comm, incoming_rank, image, add_to_front):
    mylog.debug("Receiving image from %04i", incoming_rank)
    # mylog.debug( '%04i receiving image from %04i'%(self.comm.rank,back.owner))
    arr2 = comm.recv_array(incoming_rank, incoming_rank).reshape(
        (image.shape[0], image.shape[1], image.shape[2])
    )

    if add_to_front:
        front = arr2
        back = image
    else:
        front = image
        back = arr2

    if image.shape[2] == 3:
        # Assume Projection Camera, Add
        np.add(image, front, image)
        return image

    ta = 1.0 - front[:, :, 3]
    np.maximum(ta, 0.0, ta)
    # This now does the following calculation, but in a memory
    # conservative fashion
    # image[:,:,i  ] = front[:,:,i] + ta*back[:,:,i]
    image = back.copy()
    for i in range(4):
        np.multiply(image[:, :, i], ta, image[:, :, i])
    np.add(image, front, image)
    return image


def send_to_parent(comm, outgoing_rank, image):
    mylog.debug("Sending image to %04i", outgoing_rank)
    comm.send_array(image, outgoing_rank, tag=comm.rank)


def scatter_image(comm, root, image):
    mylog.debug("Scattering from %04i", root)
    image = comm.mpi_bcast(image, root=root)
    return image
