import argparse
import json
import logging
import typing

import h5py
import gsd
import gsd.hoomd
import numpy as np

from pkg.common.args import remove_none
from pkg.common.cli import invoke_main


LOG = logging.getLogger()
GSD_SCHEMA = dict(application="", schema="hoomd", schema_version=(1, 0))
DEFAULT_BOX = (100, 100, 100)
DIMENSION = 3


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
        stage_metadata = stage_store["metadata"]
        config = json.loads(store["metadata"]["config"][()])

        with gsd.fl.open(output_filename, "w", **GSD_SCHEMA) as output:
            with gsd.hoomd.HOOMDTrajectory(output) as traj:
                match stage:
                    case "anaphase":
                        dump_trajectory(stage_store, traj, AnaphaseMod(config))
                    case "telophase":
                        dump_trajectory(stage_store, traj, TopologyMod())
                    case "relaxation" | "interphase":
                        dump_trajectory(stage_store, traj, InterphaseMod())
                    case "prometaphase":
                        dump_trajectory(stage_store, traj, PrometaphaseMod(stage_metadata))


class ParticlesData(typing.NamedTuple):
    type_ids: list[int]
    type_names: list[str]


class BondsData(typing.NamedTuple):
    pairs: list[tuple[int]]
    type_ids: list[int]
    type_names: list[str]


class TopologyMod:
    def derive_extra_particles(self, metadata: h5py.Group, next_id: int) -> ParticlesData:
        return ParticlesData([], [])

    def derive_extra_bonds(self, metadata: h5py.Group, next_id: int) -> BondsData:
        return BondsData([], [], [])

    def derive_extra_positions(self, snapshot: h5py.Group) -> np.ndarray:
        return np.zeros(shape=(0, DIMENSION))


class AnaphaseMod(TopologyMod):
    def __init__(self, config: dict):
        self._pole_position = config["mitotic_phase"]["anaphase_spindle_shift"]

    def derive_extra_particles(self, metadata: h5py.Group, next_id: int) -> ParticlesData:
        return ParticlesData(
            type_ids=[next_id],
            type_names=["spindle_pole"],
        )

    def derive_extra_bonds(self, metadata: h5py.Group, next_id: int) -> BondsData:
        pole_index = len(metadata["particle_types"])
        pairs = [(i, pole_index) for i in metadata["kinetochore_beads"]]
        return BondsData(
            pairs=pairs,
            type_ids=([next_id] * len(pairs)),
            type_names=["microtubule"],
        )

    def derive_extra_positions(self, snapshot: h5py.Group) -> np.ndarray:
        return np.reshape(self._pole_position, (1, DIMENSION))


class InterphaseMod(TopologyMod):
    def derive_extra_bonds(self, metadata: h5py.Group, next_id: int) -> BondsData:
        nucleolar_bonds = [(i, j) for i, j in metadata["nucleolar_bonds"]]
        return BondsData(
            pairs=nucleolar_bonds,
            type_ids=([next_id] * len(nucleolar_bonds)),
            type_names=["nucleolus"],
        )


class PrometaphaseMod(TopologyMod):
    def __init__(self, metadata: h5py.Group):
        self._pole_positions = metadata["pole_positions"][:]

    def derive_extra_particles(self, metadata: h5py.Group, next_id: int) -> ParticlesData:
        return ParticlesData(
            type_ids=[next_id, next_id],
            type_names=["spindle_pole"],
        )

    def derive_extra_bonds(self, metadata: h5py.Group, next_id: int) -> BondsData:
        pole_index_a = len(metadata["particle_types"])
        pole_index_b = pole_index_a + 1
        kinetochore_beads = metadata["kinetochore_beads"][:]

        pairs = []
        for chromatid_a, chromatid_b in metadata["sister_chromatids"]:
            pairs.append((kinetochore_beads[chromatid_a], pole_index_a))
            pairs.append((kinetochore_beads[chromatid_b], pole_index_b))

        return BondsData(
            pairs=pairs,
            type_ids=([next_id] * len(pairs)),
            type_names=["microtubule"],
        )

    def derive_extra_positions(self, snapshot: h5py.Group) -> np.ndarray:
        return self._pole_positions


def dump_trajectory(
    store: h5py.Group,
    traj: gsd.hoomd.HOOMDTrajectory,
    topology_mod: TopologyMod,
):
    metadata = store["metadata"]
    particles = derive_particles(metadata, topology_mod)
    bonds = derive_bonds(metadata, topology_mod)
    box_shape = DEFAULT_BOX

    # Dump all snapshots.
    for step in store[".steps"]:
        snapshot = store[step]

        stored_positions = snapshot["positions"][:]
        extra_positions = topology_mod.derive_extra_positions(snapshot)
        positions = np.concatenate([stored_positions, extra_positions])

        frame_data = FrameData(
            step=int(step),
            box_shape=box_shape,
            particle_types=particles.type_ids,
            particle_type_names=particles.type_names,
            particle_positions=positions,
            particle_attributes={},
            bond_types=bonds.type_ids,
            bond_type_names=bonds.type_names,
            bond_pairs=bonds.pairs,
        )
        traj.append(make_hoomd_frame(frame_data))


def derive_particles(metadata: h5py.Group, topology_mod: TopologyMod) -> ParticlesData:
    stored_types = metadata["particle_types"][:]
    stored_type_names = [
        name for name, _type_id in sorted(
            stored_types.dtype.metadata["enum"].items(),
            key=(lambda entry: entry[1]),
        )
    ]

    extra_particles = topology_mod.derive_extra_particles(metadata, next_id=len(stored_type_names))

    return ParticlesData(
        type_ids=(list(stored_types) + extra_particles.type_ids),
        type_names=(stored_type_names + extra_particles.type_names),
    )


def derive_bonds(metadata: h5py.Group, topology_mod: TopologyMod) -> BondsData:
    chain_ranges = metadata["chain_ranges"][:]

    stored_pairs = sum(
        (define_linear_bonds(start, end) for start, end in chain_ranges),
        [],
    )
    stored_type_ids = [0] * len(stored_pairs)
    stored_type_names = ["chrom"]

    extra_bonds = topology_mod.derive_extra_bonds(metadata, next_id=len(stored_type_names))

    return BondsData(
        pairs=(stored_pairs + extra_bonds.pairs),
        type_ids=(stored_type_ids + extra_bonds.type_ids),
        type_names=(stored_type_names + extra_bonds.type_names),
    )


def define_linear_bonds(start: int, end: int) -> list[tuple[int, int]]:
    return list(zip(range(start, end - 1), range(start + 1, end)))


class FrameData(typing.NamedTuple):
    step: int
    box_shape: tuple[float, float, float]
    particle_type_names: list[str]
    particle_types: np.ndarray
    particle_positions: np.ndarray
    particle_attributes: dict[str, np.ndarray]
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

    # Custom attributes
    for key, values in data.particle_attributes.items():
        frame.log[f"particles/{key}"] = values

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
