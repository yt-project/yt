import numpy as np

from .data_structures import GadgetBinaryHeader, GadgetDataset
from .definitions import gadget_field_specs, gadget_ptype_specs
from .io import IOHandlerGadgetBinary

vector_fields = dict(IOHandlerGadgetBinary._vector_fields)

block_ids = {
    "Coordinates": "POS",
    "Velocities": "VEL",
    "ParticleIDs": "ID",
    "Mass": "MASS",
    "InternalEnergy": "U",
    "Density": "RHO",
    "SmoothingLength": "HSML",
}


def write_record(fp, data, endian):
    dtype = endian + "i4"
    size = np.array(data.nbytes, dtype=dtype)
    fp.write(size.tobytes())
    fp.write(data.tobytes())
    fp.write(size.tobytes())


def write_block(fp, data, endian, fmt, block_id):
    assert fmt in [1, 2]
    block_id = "%-4s" % block_id
    if fmt == 2:
        block_id_dtype = np.dtype([("id", "S", 4), ("offset", endian + "i4")])
        block_id_data = np.zeros(1, dtype=block_id_dtype)
        block_id_data["id"] = block_id
        block_id_data["offset"] = data.nbytes + 8
        write_record(fp, block_id_data, endian)
    write_record(fp, data, endian)


def fake_gadget_binary(
    filename="fake_gadget_binary",
    npart=(100, 100, 100, 0, 100, 0),
    header_spec="default",
    field_spec="default",
    ptype_spec="default",
    endian="",
    fmt=2,
):
    """Generate a fake Gadget binary snapshot."""
    header = GadgetBinaryHeader(filename, header_spec)
    field_spec = GadgetDataset._setup_binary_spec(field_spec, gadget_field_specs)
    ptype_spec = GadgetDataset._setup_binary_spec(ptype_spec, gadget_ptype_specs)
    with open(filename, "wb") as fp:
        # Generate and write header blocks
        for i_header, header_spec in enumerate(header.spec):

            specs = []
            for name, dim, dtype in header_spec:
                # workaround a FutureWarning in numpy where np.dtype(name, type, 1)
                # will change meaning in a future version so
                name_dtype = [name, endian + dtype, dim]
                if dim == 1:
                    name_dtype.pop()
                specs.append(tuple(name_dtype))

            header_dtype = np.dtype(specs)
            header = np.zeros(1, dtype=header_dtype)
            if i_header == 0:
                header["Npart"] = npart
                header["Nall"] = npart
                header["NumFiles"] = 1
                header["BoxSize"] = 1
                header["HubbleParam"] = 1
            write_block(fp, header, endian, fmt, "HEAD")

        npart = dict(zip(ptype_spec, npart))
        for fs in field_spec:
            # Parse field name and particle type
            if isinstance(fs, str):
                field = fs
                ptype = ptype_spec
            else:
                field, ptype = fs
                if isinstance(ptype, str):
                    ptype = (ptype,)
            # Determine field dimension
            if field in vector_fields:
                dim = vector_fields[field]
            else:
                dim = 1
            # Determine dtype (in numpy convention)
            if field == "ParticleIDs":
                dtype = "u4"
            else:
                dtype = "f4"
            dtype = endian + dtype
            # Generate and write field block
            data = []
            for pt in ptype:
                data += [np.random.rand(npart[pt], dim)]
            data = np.concatenate(data).astype(dtype)
            if field in block_ids:
                block_id = block_ids[field]
            else:
                block_id = ""
            write_block(fp, data, endian, fmt, block_id)
    return filename
