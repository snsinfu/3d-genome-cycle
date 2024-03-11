# Simulation parameters


## config.json

The driver script (scripts/simulate) takes a JSON file (config.json) for
configuring technical parameters used in the simulation. The JSON file
should contain two subobjects, "anatelophase" and "interphase":

```json
{
    "anatelophase": {
        "coarse_graining": 100,
        ...
    },
    "interphase": {
        "temperature": 1,
        ...
    }
}
```

The tables below summarize the configuration keys in each subobject. All keys
are optional and the indicated default values are used if not specified.


### anatelophase

| key                 | default   | description |
|---------------------|-----------|-------------|
| `coarse_graining`   | 100       | The level of further coarse-graining in anatelophase simulation. |
| `temperature`       | 1         | The temperature (kT) for the overdamped Langevin dynamics. |
| `timestep`          | 1e-4      | The size of discretized timestep in hours. |
| `dragging_steps`    | 200000    | The number of timesteps in the dragging stage. |
| `packing_steps`     | 50000     | The number of timesteps in the packing stage. |
| `sampling_interval` | 1000      | The interval, in timesteps, of which snapshots are taken and saved to the output file. |
| `logging_interval`  | 10000     | The interval, in timesteps, of which the progress of the simulation is printed. |
| `start_center`      | [0, 0, 5] | The coordinates, in μm, of the point where chromosomes are scattered around before the dragging simulation starts. |
| `start_stddev`      | 1         | The standard deviation of the distribution of chromosomes around start_center. |
| `core_diameter`     | 0.3       | The diameter, in μm, of beads that represent coarse-grained chromosome segments. |
| `core_repulsion`    | 2         | The repulsive energy of beads. |
| `bond_length`       | 0.3       | The bond length of the beads in μm. |
| `bond_spring`       | 1000      | The spring constant in kT/μm². |
| `bending_energy`    | 0         | The energetic cost of bending, in kT, applied to the chromosomes. |
| `core_mobility`     | 1         | The mobility of beads in μm²/(kT×hour). |
| `dragging_spring`   | 1.5       | The spring constant of the dragging force applied to the kinetochore of the chromosomes. |
| `packing_radius`    | 1.2       | The radius of the packing sphere, in μm, used in the packing stage. |
| `packing_spring`    | 10        | The spring constant to drag chromosomes into the packing sphere, in kT/μm². |


### interphase

| key                          | default         | description |
|------------------------------|-----------------|-------------|
| `temperature`                | 1               | The tempearture (kT) for the overdamped Langevin dynamics. |
| `timestep`                   | 1e-5            | The size of discretized timestep in hours. |
| `steps`                      | 700000          | The number of timesteps in the interphase simulation |
| `relaxation_steps`           | 10000           | The number of timesteps in the initial relaxation. |
| `relaxation_spacestep`       | 0.001           | The maximum movement of beads within a single relaxation timestep in μm, to avoid numerical instability in the initial, high-energy conformation. |
| `sampling_interval`          | 1000            | The interval, in timesteps, of which snapshots are taken and saved to the output file. |
| `logging_interval`           | 10000           | The interval, in timesteps, of which the progress of the simulation is printed. |
| `contactmap_distance`        |
| `contactmap_update_interval` | 20              | The interval, in timesteps, of which instantaneous contact map is calculated. |
| `contactmap_output_window`   | 10              | The interval, in frames (sampling_interval), of which the time-integrated contact map is saved to the output file. |
| `a_core_diameter`            | 0.30            | The diameter of type-A beads in μm. |
| `b_core_diameter`            | 0.24            | The diameter of type-B beads in μm. |
| `a_core_repulsion`           | 2.5             | The repulsive energy of type-A beads in kT. |
| `b_core_repulsion`           | 2,5             | The repulsive energy of type-B beads in kT. |
| `a_core_bond_spring`         | 100             | The bond spring constant for type-A beads in kT/μm². |
| `b_core_bond_spring`         | 50              | The bond spring constant for type-B beads in kT/μm². |
| `a_core_bond_length`         | 0               | The bond length between adjacent type-A beads on the chain. |
| `b_core_bond_length`         | 0               | The bond length between adjacent type-B beads on the chain. |
| `a_core_mobility`            | 1               | The mobility of type-A beads in μm²/(kT hour). |
| `b_core_mobility`            | 1               | The mobility of type-B beads in μm²/(kT hour). |
| `core_scale_init`            | 0.5             | The initial scaling factor for the bead diameter. |
| `core_scale_tau`             | 0.5             | The time constant, in hour, for the bead-diameter scaling. |
| `bond_scale_init`            | 0.5             | The initial scaling factor fot the bond length. |
| `bond_scale_tau`             | 0.5             | The time constant, in hour, for the bond-length scaling. |
| `nucleolus_bead_count`       | 2               | The number of virtual nucleolar particles that bind to a chromatin bead in NOR. |
| `nucleolus_ab_factor`        | [0, 10]         | The (A, B) repulsive interaction weights for the virtual nucleolar particles. |
| `nucleolus_bond_spring`      | 10              | The spring constant for the bond between chromatin beads and nucleolar particles. |
| `nucleulus_bond_length`      | 0               | The length for the bond between chromatin beads and nucleolar particles. |
| `nucleolus_droplet_energy`   | 0.3             | Attractive interaction energy, in kT, among nucleolar particles. |
| `nucleolus_droplet_decay`    | 0.2             | Distance scale, in μm, of the decaying tail of the attractive interaction among nucleolar particles. |
| `nucleolus_droplet_cutoff`   | 0.4             | The cutoff distance, in μm, of the attractive interaction among nucleolar particles. |
| `nucleolus_mobility`         | 1               | The mobility of nucleolar particles in μm²/(kT hour). |
| `wall_semiaxes_init`         | [2, 2, 2]       | The initial lengths, in μm, of the three axes of the ellipsoidal wall. |
| `wall_semiaxes_spring`       | [3e4, 3e4, 3e4] | The spring constants applied to the three axes in μm²/(kT hour). |
| `wall_packing_spring`        | 1000            | The spring constant for constraining all particles within the wall. |
| `wall_ab_factor`             | [0, 10]         | The (A, B) repulsive interaction weights for the internal side of the wall. |
| `wall_mobility`              | 2e-4            | The mobility of the whole wall in μm²/(kT hour). |

