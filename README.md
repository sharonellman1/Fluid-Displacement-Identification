# Fluid-Displacement-Identification

This is fluid displacement identification code adapted from the code used in Mascini et al., 2020. It can be used to determine properties for individual displacements, such as time at which a displacement occurred, which pores were involved, contact angle per displacement, etc. When working with imbibition datasets, it also gives information over event I_n type and snap-offs. 

It requires a pore network extraction in statoil format. The maximum inscribed spheres of the pore network extraction are then overlayed on the segmented fluid distribution images to determine fluid occupancy, contact angle and curvature per pore/throat. As input it requires:
- link and node files in statoil format
- majority and mean (oil = 1, brine+rock = 0) in pore maximum inscribed spheres per timestep
- majority and mean in throat maximum inscribed spheres per timestep
- majority and mean on throat surfaces per timestep
- the volumes of each pore maximum inscribed sphere
- Pc over the entire field of view per timestep
- average curvature over the entire field of view per timestep
- curvature and contact angle per pore per timestep
- curvature and contact angle over the entire field of view in the final timestep

The file "networkComponents.m" used in "eventClusteringAnalysis.m" can be downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/42040-find-network-components?s_tid=srchtitle

## References
- Mascini, A., Cnudde, V. and Bultreys, T., 2020. Event-based contact angle measurements inside porous media using time-resolved micro-computed tomography. Journal of colloid and interface science, 572, pp.354-363. https://doi.org/10.1016/j.jcis.2020.03.099
- Ellman, S., Mascini, A., Bultreys, T., 2023. Validating mechanistic models of fluid displacement during imbibition. Advances in Water Resources. https://doi.org/10.1016/j.advwatres.2023.104590
