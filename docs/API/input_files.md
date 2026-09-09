---
icon: lucide/file-input
---

# Input files

## JSON input file format

Many CabanaPD numerical inputs are specified via JSON values.
Others, such as the model to be used can only be specified in `.cpp` files.

All CabanaPD examples are within the `examples/` directory; subdirectories contain the
primary categories: `dem`, `mechanics`, and `thermomechanics`, and each of these contains
an `inputs/` subdirectory with the JSON input files.

### Required inputs

The following are required for any CabanaPD simulation:

 - System size `system_size` or both `low_corner` and `high_corner` (per dimension)
 - Number of cells `num_cells` or spacing `dx` (per dimension)
 - Mass density `density`
 - Bulk modulus `bulk_modulus` or elastic modulus `elastic_modulus` (potentially per material)
 - Neighbor horizon (search cutoff) `horizon` or m-ratio `m` (horzizon / dx)
 - Timestep `timestep`
 - Final simulation time `final_time`
 - Output step frequency `output_frequency`

Further inputs are commonly used across example cases, but not required.

## Mapping input values to variables in the simulation

To extract values from JSON inputs within a given case, use:

```
CabanaPD::Inputs inputs( json_file );
auto variable = inputs["variable"];
```
