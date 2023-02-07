# TO-DO
Upgrades:
* Upgrade from Julia 1.8.2 to Julia 1.8.5
* Upgrade all packages to last version.

Dashboard/launch & Worker:
* Add config option for launching workers with a built sysimage
* Add a way of launching and visualizing individual simulations.

TumorModel:
* Re-check how each of the properties (pr/dr...) is modelled to make sure they are OK.

Simulate:
* Add a column for "n0 % to carrying capacity" -> t_detecting_size / s_size - In other papers its mentioned as initial size % to carrying capacity.

Analysis:
* Add a function to calculate from genotype, treatment and t_detecting_size the % of resistant cells on detection.

Results:
* Check if the results so far match the ones from the bibliography. In particular  (A.R. Strobl, 2021) and (A.R. Strobl, 2022)