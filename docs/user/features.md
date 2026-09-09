---
icon: lucide/rocket
---

# CabanaPD features

CabanaPD currently includes the following:

  - Force models
    - PD bond-based (pairwise): PMB (prototype microelastic brittle)
    - PD state-based (many-body): LPS (linear peridynamic solid)
    - DEM (contact): normal repulsion, Hertzian, HertzianJKR (Johnson–Kendall–Roberts)
    - Hybrid PD-DEM
    - Multi-material systems can be constructed for any models of the **same category**
      (bond-based, state-based, contact) above. Cross-term interactions can be averaged, requiring **identical model** types
  - Mechanical response
    - Elastic only (no failure)
    - Brittle fracture
    - Elastic-perfectly plastic (*Currently bond-based only*)
  - Thermomechanics (*Currently bond-based only, single material only*)
    - Optional heat transfer
  - Time integration
    - Velocity Verlet
  - Pre-crack creation
  - Particle boundary conditions
    - Body terms which apply to all particles
  - Grid-based particle generation supporting custom geometry
  - Output options
    - Total strain energy density
    - Total damage (if fracture is enabled)
    - Per particle output using HDF5 or SILO
        - Base fields: position (reference or current), displacement, velocity, force, material type
        - Strain energy density, damage
        - Stress
        - LPS fields (if used): weighted volume, dilatation
        - Thermal fields (if used): temperature
