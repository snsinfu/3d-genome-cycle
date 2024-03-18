import argparse
import dataclasses
import logging

import h5py
import gsd
import gsd.hoomd
import numpy as np

from pkg.common.args import remove_none
from pkg.common.cli import invoke_main


LOG = logging.getLogger()
GSD_SCHEMA = dict(application="", schema="hoomd", schema_version=(1, 0))


def main(
    *,
    input_filename: str,
    output_filename: str,
    stage: str = "interphase",
):
    # HOOMDTrajectory logs confusing messages at the INFO level.
    logging.getLogger("gsd.hoomd").setLevel(logging.WARNING)

    with h5py.File(input_filename, "r") as store:
        stage_store = store["stages"][stage]

        with gsd.fl.open(output_filename, "w", **GSD_SCHEMA) as output:
            with gsd.hoomd.HOOMDTrajectory(output) as traj:
                match stage:
                    case "anaphase" | "telophase":
                        dump_metaphase(stage_store, traj)
                    case "relaxation" | "interphase":
                        dump_interphase(stage_store, traj)


def dump_metaphase(
    store: h5py.Group,
    traj: gsd.hoomd.HOOMDTrajectory,
):
    metadata = store["metadata"]
    particle_types = metadata["particle_types"][:]
    chain_ranges = metadata["chain_ranges"][:]
    kinetochore_beads = metadata["kinetochore_beads"][:]

    # GSD expects the type IDs to be the indices into this list.
    particle_type_names = [
        name for name, _type_id in sorted(
            particle_types.dtype.metadata["enum"].items(),
            key=(lambda entry: entry[1]),
        )
    ]

    # Add a virtual spindle pole body particle at the origin for visualizing
    # the dragging force (or microtubules).
    spb_position = (0, 0, 0)
    spb_index = len(particle_types)
    spb_type_id = len(particle_type_names)
    particle_types = np.concatenate([particle_types, [spb_type_id]])
    particle_type_names.append("spb")

    spb_bonds = [(kin_index, spb_index) for kin_index in kinetochore_beads]

    # Deduce implicit bonds.
    chromosomal_bonds = np.concatenate([
        define_linear_bonds(start, end)
        for start, end in chain_ranges
    ])
    bond_pairs = np.concatenate([
        chromosomal_bonds,
        spb_bonds,
    ])
    bond_types = np.concatenate([
        [0] * len(chromosomal_bonds),
        [1] * len(spb_bonds),
    ])
    bond_type_names = ["chrom", "spb"]

    # Our system is open, so define an arbitrary box. Note that OVITO assumes
    # that the system is periodic and wraps bonds at the boundaries, so the box
    # needs to be large.
    box_shape = (100, 100, 100)

    # Dump all snapshots.
    for step in store[".steps"]:
        snapshot = store[step]
        positions = np.concatenate([snapshot["positions"][:], [spb_position]])

        frame_data = FrameData(
            step=int(step),
            box_shape=box_shape,
            particle_type_names=particle_type_names,
            particle_types=particle_types,
            particle_positions=positions,
            bond_type_names=bond_type_names,
            bond_pairs=bond_pairs,
            bond_types=bond_types,
        )
        traj.append(make_hoomd_frame(frame_data))


def dump_interphase(
    store: h5py.Group,
    traj: gsd.hoomd.HOOMDTrajectory,
):
    metadata = store["metadata"]
    particle_types = metadata["particle_types"][:]
    chain_ranges = metadata["chain_ranges"][:]
    nucleolar_bonds = metadata["nucleolar_bonds"][:]

    # GSD expects that the type IDs are the indices into this list.
    particle_type_names = [
        name for name, _type_id in sorted(
            particle_types.dtype.metadata["enum"].items(),
            key=(lambda entry: entry[1]),
        )
    ]

    # Deduce implicit bonds.
    chromosomal_bonds = np.concatenate([
        define_linear_bonds(start, end)
        for start, end in chain_ranges
    ])
    bond_pairs = np.concatenate([
        chromosomal_bonds,
        nucleolar_bonds,
    ])
    bond_types = np.concatenate([
        [0] * len(chromosomal_bonds),
        [1] * len(nucleolar_bonds),
    ])
    bond_type_names = ["chrom", "nucleo"]

    # Our system is open, so define an arbitrary box. Note that OVITO assumes
    # that the system is periodic and wraps bonds at the boundaries, so the box
    # needs to be large.
    # - Should we estimate the enclosing box?
    # - How can we translate the center of the box to the origin?
    box_shape = (100, 100, 100)

    # Dump all snapshots.
    for step in store[".steps"]:
        snapshot = store[step]
        positions = snapshot["positions"][:]

        frame_data = FrameData(
            step=int(step),
            box_shape=box_shape,
            particle_type_names=particle_type_names,
            particle_types=particle_types,
            particle_positions=positions,
            bond_type_names=bond_type_names,
            bond_pairs=bond_pairs,
            bond_types=bond_types,
        )
        traj.append(make_hoomd_frame(frame_data))


def define_linear_bonds(start: int, end: int) -> np.ndarray:
    indices = np.arange(start, end)
    return np.transpose([indices[:-1], indices[1:]])


@dataclasses.dataclass
class FrameData:
    step: int
    box_shape: tuple[float, float, float]
    particle_type_names: list[str]
    particle_types: np.ndarray
    particle_positions: np.ndarray
    bond_type_names: list[str]
    bond_pairs: np.ndarray
    bond_types: np.ndarray


def make_hoomd_frame(data: FrameData) -> gsd.hoomd.Frame:
    frame = gsd.hoomd.Frame()

    frame.configuration = gsd.hoomd.ConfigurationData()
    frame.configuration.box = (*data.box_shape, 0, 0, 0)
    frame.configuration.step = data.step

    frame.particles = gsd.hoomd.ParticleData()
    frame.particles.types = data.particle_type_names
    frame.particles.position = data.particle_positions
    frame.particles.typeid = data.particle_types
    frame.particles.N = len(frame.particles.position)

    frame.bonds = gsd.hoomd.BondData(2)
    frame.bonds.types = data.bond_type_names
    frame.bonds.group = data.bond_pairs
    frame.bonds.typeid = data.bond_types
    frame.bonds.N = len(frame.bonds.group)

    # These entries are unused but need to be set.
    # See: https://gsd.readthedocs.io/en/v3.2.1/python-module-gsd.hoomd.html#gsd.hoomd.BondData
    frame.angles = gsd.hoomd.BondData(3)
    frame.dihedrals = gsd.hoomd.BondData(4)
    frame.impropers = gsd.hoomd.BondData(4)
    frame.pairs = gsd.hoomd.BondData(2)
    frame.constraints = gsd.hoomd.ConstraintData()

    frame.state = {}
    frame.log = {}

    frame.validate()

    return frame


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=str)
    parser.add_argument("input_filename", type=str)
    parser.add_argument("output_filename", type=str)
    return remove_none(vars(parser.parse_args()))


if __name__ == "__main__":
    invoke_main(main, parse_args(), LOG)
