---
title: A Sea-Breeze Parametrisation
permalink: /
---
This is a  documentation of a parametrization of sub-grid scale lan-sea-breeze 
circulation systems. The here described method can be applied in numerical 
weather prediction and climate model to inform cumulus parametrisations about 
the presence of any sub-grid scale sea-breezes or just to calculated global 
climatologies of sea-breeze occurances.
## Background
The overall dependence of the formation of clouds and rainfall near coasts on large-scale weather patterns is weaker than over the open ocean or inland areas
    [Bergemann and Jakob (2016)](https://arxiv.org/abs/1603.02392v1).
    This is due to the fact that clouds and rainfall are often influenced by land-sea-breezes; 
    local winds that blow on and off-shore along the coastline. As these breezes 
    are not included in current weather and climate models rainfall near coasts 
    is poorly simulated in these models. This GitHub page presents an algorithm, 
    or trigger, that is used as a parametrisation fo land-sea-breezes. 
    The presented code can be applied to inform a cumulus parametrization scheme 
    model to increase the occurrence of convective clouds when land-see-breezes 
    are diagnosed to be important. A detailed scientific description of the 
    algorithm can be found in [Bergemann et al. (2017)](http://onlinelibrary.wiley.com/doi/10.1002/2017MS001048/full)
   

### The Algorithm
<p>To inform the global model about the presence of meso-scale systems
a simple filter algorithm is applied. The filter is based on large-scale
atmospheric condtions in a climate model (g) and consideres six criteria (see also Fig.1):
<ol><li>Change of large-scale wind direction &Delta;&alpha;</li>
<li>Change in large-scale wind speed |&Delta;V|</li>
<li>Magnitude of large-scale wind |V|</li>
<li>Temperature contrast between land and ocean |&Delta;T|</li></ol></p>
<figure>
<img src="assets/images/Seabreeze_detect.png">
<figcaption>Fig.1 - Schematic of the filtering process to identify potential sea-breeze conditions. 
The blue letters represent applied thresholds, the g's indicate the values of the considered variable from the model.
</figcaption>
</figure>
### Filtering of Coastal Areas
The above described algorithm is only applied in coastal areas that are
influenced by land-sea-breeze circulation systems. This radius of influence
is typically in the order of &cong; 150km on- and offshore along the coastline.
The coastline is calculated by applying a [_Sobel Operator_](https://en.wikipedia.org/wiki/Sobel_operator) 
to the models land area fraction data. The sea-ice 
fraction is also taken into account when calculating the coastline. 
This is done to avoid the application of the sea-breeze algorithm over areas 
covered with sea-ice.
### Steep Terrain
Coastal areas with steep mountain terrain can cause problems to the algorithm,
especially on a coars resolution (&ge 75km). In areas with steep topography, 
large standard deviation of the sub-grid scale orogaphy field, a moist adiabatic
decent to sea-level hieght (z<sub>0</sub>) is performed and the theoretical 
temperature at z<sub>0</sub> is calculated. To keep the temperature field 
continous the sub-grid orography field is filtered by a [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function) and a simple moist adiabitic 
decent is performed to the field:

<center>T = T<sub>z</sub> - ( &Gamma;<sub>moist</sub> &sdot; z &sdot; sigmoid(&sigma;<sub>oro</sub>) )</center>

## Implementation
The parametrisation is implemented using the fortran 90 standard. The source code file :
```
sea_breeze_diag.f90
```
which contains the module *sea_breeze_diag_mod* this module has the following subroutines:
* seabreeze_diag:
  - Contains the sea-breeze filtering process
* get_edges:
  - Applies the Sobel operator to create a coast line
* get_dist:
  - Calculates the Euclidian distance to the closest coastline point calculated by *get_edges*
*get_dist* is only called by *get_edges*.

All subroutines make use of information from neighboring grid points. It is therefore important
that an inter processor communication for boundary swapping (e.g swap_bounds in case of UM) is
established before calling the subroutines. 

If you are considering implementing this routine take a look at the  [about](/about) 
section to learn more on how to contribute and improve this project. You are encouraged 
to get in touch via [GitHub](https://github.com/antarcticrainforest/seabreeze_param). 
Bugs should be reported either on the GitHub [issues](https://github.com/antarcticrainforest/seabreeze_param/issues) 
pages or by sending an [email](mailto:martin.bergemann@monash.edu) to the author of this page.
