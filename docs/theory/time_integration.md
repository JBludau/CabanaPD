---
icon: lucide/clock
---

# Integration of the particle displacement in time

To integrate the [equation of motion](peridynamics.md#theoretical-framework) in time CabanaPD provides two integration schemata.
The [velocity verlet](time_integration.md#velocity-verlet-for-non-stationary-evolution-in-time) is an explicit split step integration in time.
It should be used whenever CabanaPD simulates non-stationary processes.
When searching for stationary points, the adaptive dynamic relaxation (ADR) scheme can speed up the computation significantly.
Crucially, both time integration schemes can be used interchangeably.
An example that uses ADR for the quasi-stationary part of the simulation and velocity verlet for the failure is the [`dogboneADR`](../user/examples.md#mechanics).

## Velocity-Verlet for non-stationary evolution in time

## Adaptive Dynamic Relaxation for finding stationary solutions


## References
