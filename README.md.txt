This is fluid displacement identification code adapted from the code used in Mascini et al., 2020. It can be used to determine properties for individual displacements, such as time at which a displacement occurred, which pores were involved, contact angle per displacement, etc. When working with imbibition datasets, it also gives information over event I_n type and snap-offs. 

It requires a pore network extraction in statoil format. The maximum inscribed spheres of the pore network extraction are then overlayed on the segmented fluid distribution images to determine fluid occupancy, contact angle and curvature per pore/throat. As input it requires:
- link and node files in statoil format
- majority and mean (oil = 1, brine+rock = 0) in pore maximum inscribed spheres
- majority and mean in throat maximum inscribed spheres
- majority and mean on throat surfaces
- the volumes of each pore maximum inscribed sphere
- Pc over the entire field of view per timestep
- average curvature over the entire field of view per timestep
- curvature and contact angle per pore per timestep
- curvature and contact angle over the entire field of view in the final timestep



 