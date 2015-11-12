# TaSMET
Thermoacoustic System Modeling Environment Twente

Welcome to the ThermoAcoustic System Modeling Environment Twente, or
TASMET. TASMET is a computer code to model thermoacoustic (TA)
engines, refrigerators and combined systems by providing nonlinear
models for laminar/turbulent oscillating flow in ducts, heat
exchangers and stacks/regenerators. A coupling to the mechanical and
electrical domain is provided with a piston model.

- For a quick glance on the possibilities of TaSMET, please look into
the User's guide provided with the code.
- For the installation guide, please look into the [Wiki][Wiki]
- Examples of working models will be provided. We are currently
  working on that. Once a regenerator model is implemented, we will
  provide examples of Stirling models.


[Wiki]: https://github.com/amdj/tasmet/wiki

## Example models

Example models for TaSMET will be provided in the
[tasmet_examples][ta_ex] repository.



[ta_ex]: https://github.com/amdj/tasmet_examples


## Release notes version 0.1
- Stacks cannot yet be connected to ConnectorVolumes
- Stacks cannot yet be connected to Pistons
- VelocityBc probably still has a bug

