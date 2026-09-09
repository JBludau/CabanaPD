---
icon: lucide/book-open-check
---

## Mechanics
Examples which only include mechanics and fracture are within `examples/mechanics`.

 -  The first example is an elastic wave propagating through a cube from an initial Gaussian radial displacement profile[^1]. The example can be run with:
    ```
    ./CabanaPD/build/install/bin/ElasticWave CabanaPD/examples/mechanics/inputs/elastic_wave.json
    ```

 -  The next example is the Kalthoff-Winkler experiment[^2], where an impactor causes crack propagation at an angle from two pre-notches on a steel plate.
    ```
    ./CabanaPD/build/install/bin/KalthoffWinkler CabanaPD/examples/mechanics/inputs/kalthoff_winkler.json
    ```

 -  Another example is crack branching in a pre-notched soda-lime glass plate due to traction loading[^3].
    ```
    ./CabanaPD/build/install/bin/CrackBranching CabanaPD/examples/mechanics/inputs/crack_branching.json
    ```

 -  An example with multiple random pre-notches is also available.
    ```
    ./CabanaPD/build/install/bin/RandomCracks CabanaPD/examples/mechanics/inputs/random_cracks.json
    ```

 -  The next example is a fragmenting cylinder due to internal pressure[^4].
    This problem can either run with PD only or with hybrid PD-DEM contact.
    ```
    ./CabanaPD/build/install/bin/FragmentingCylinder CabanaPD/examples/mechanics/inputs/fragmenting_cylinder.json
    ```

 -  An example highlighting plasticity simulates a tensile test based on an ASTM standard dogbone specimen.
    ```
    ./CabanaPD/build/install/bin/DogboneTensileTest CabanaPD/examples/mechanics/inputs/dogbone_tensile_test.json
    ```

 -  An example demonstrating the peridynamic stress tensor computation simulates a square plate under tension with a circular hole at its center[^5].
    ```
    ./CabanaPD/build/install/bin/PlateWithHole CabanaPD/examples/mechanics/inputs/plate_with_hole.json
    ```

 -  An example of multi-material simulation demonstrates crack propagation in a pre-notched plate with a stiff inclusion under traction loading.
    ```
    ./CabanaPD/build/install/bin/CrackInclusion CabanaPD/examples/mechanics/inputs/crack_inclusion.json
    ```

## Powder dynamics
Examples which only include mechanics and fracture are within `examples/dem`.

 -  An example using DEM-only demonstrates powder filling in a container.
    ```
    ./CabanaPD/build/install/bin/PowderFill CabanaPD/examples/mechanics/inputs/powder_fill.json
    ```

## Thermomechanics
Examples which demonstrate temperature-dependent mechanics and fracture are within `examples/thermomechanics`.

 -  The first example is thermoelastic deformation in a homogeneous plate due to linear thermal loading[^6].
    ```
    ./CabanaPD/build/install/bin/ThermalDeformation CabanaPD/examples/thermomechanics/thermal_deformation.json
    ```

 -  The second example is crack initiation and propagation in an alumina ceramic plate due to a thermal shock caused by water quenching[^7].
    ```
    ./CabanaPD/build/install/bin/ThermalCrack CabanaPD/examples/thermomechanics/thermal_crack.json
    ```

## Thermomechanics with heat transfer
Examples with heat transfer are within `examples/thermomechanics`.

 -  The first example is pseudo-1d heat transfer (no mechanics) in a cube.
    ```
    ./CabanaPD/build/install/bin/ThermalDeformationHeatTransfer CabanaPD/examples/thermomechanics/heat_transfer.json
    ```
    The same example with fully coupled thermomechanics can be run (with a much smaller timestep) using `thermal_deformation_heat_transfer.json`.

 -  The second example is pseudo-1d heat transfer (no mechanics) in a pre-notched cube.
    ```
    ./CabanaPD/build/install/bin/ThermalDeformationHeatTransferPrenotched CabanaPD/examples/thermomechanics/heat_transfer.json
    ```


## References

[^1]: P. Seleson and D.J. Littlewood, Numerical tools for improved convergence of meshfree peridynamic discretizations, in Handbook of Nonlocal Continuum Mechanics for Materials and Structures, G. Voyiadjis, ed., Springer, Cham, 2018. [doi:10.1007/978-3-319-22977-5_39-1](https://doi.org/10.1007/978-3-319-22977-5_39-1)

[^2]: J.F. Kalthoff and S. Winkler, Failure mode transition at high rates of shear loading, in Impact Loading and Dynamic Behavior of Materials, C.Y. Chiem, H.-D.
    Kunze, and L.W. Meyer, eds., Vol 1, DGM Informationsgesellschaft Verlag (1988) 185-195.

[^3]: F. Bobaru and G. Zhang, Why do cracks branch? A peridynamic investigation of dynamic brittle fracture, *International Journal of Fracture* **196** (2015): 59–98. [doi:10.1007/s10704-015-0056-8](https://doi.org/10.1007/s10704-015-0056-8)

[^4]: D.J. Littlewood, M.L. Parks, J.T. Foster, J.A. Mitchell, and P. Diehl, The peridigm meshfree peridynamics code, *Journal of Peridynamics and Nonlocal Modeling* **6** (2024): 118–148. [doi:10.1007/s42102-023-00100-0](https://doi.org/10.1007/s42102-023-00100-0)

[^5]: A.S. Fallah, I.N. Giannakeas, R. Mella, M.R. Wenman, Y. Safa, and H. Bahai, On the computational derivation of bond-based peridynamic stress tensor, *Journal of Peridynamics and Nonlocal Modeling* **2** (2020): 352–378. [doi:10.1007/s42102-020-00036-9](https://doi.org/10.1007/s42102-020-00036-9)

[^6]: D. He, D. Huang, and D. Jiang, Modeling and studies of fracture in functionally graded materials under thermal shock loading using peridynamics, *Theoretical and Applied Fracture Mechanics* **111** (2021): 102852. [doi:10.1016/j.tafmec.2020.102852](https://doi.org/10.1016/j.tafmec.2020.102852)

[^7]: C.P. Jiang, X.F. Wu, J. Li, F. Song, Y.F. Shao, X.H. Xu, and P. Yan, A study of the mechanism of formation and numerical simulations of crack patterns in ceramics subjected to thermal shock, *Acta Materialia* **60**(11) (2012): 4540–4550. [doi:10.1016/j.actamat.2012.05.020](https://doi.org/10.1016/j.actamat.2012.05.020)
