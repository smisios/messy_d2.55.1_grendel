Examples for clams-2.0-messy driven by era-interim data

* 2.0_era
CLaMS-1.0 configuartion for climat. run

* 2.0_era_l1.5
CLaMS-2.0 configuartion for climat. run with:
Delta t=24h, l=1.5 1/day, fac_eliminate=1.0 - like in CLaMS-1.0
Thickness of the lower layer: 70K (instead of 50K)
new init-files (r=100km in the lowest layer for all resolutions)
TR/DC tracer instead of TRACER

* 2.0_era_l4.0_nodc 
Delta t=6h, l=4.0 1/day, fac_eliminate=0.3

* 2.0_era_l4.0_dc
deep convection is switched on with:
cbvf=0.0 (1/s^2),   (critical wet BVF)
crate=1.,   (no overshooting)
min_dtheta=10., (minimal Delta theta for switching on deep conv)

* 2.0_era_l4.0_d0.0/_perp
Vertical mixing is included with dry BVF=0.0 (1/s^2)
(perp) configuration for perp-run

* 2.0_era_l4.0_d0.5
Vertical mixing is enhanced to BVF=0.5e-4 (1/s^2)

Copy nml-lists from every dir to out_1/out_2/out_fast/out_perp
and adjust respective files in:

/p/home/jusers/konopka1/juwels/jicg1119/messy-clams/messy/util

(xmessy_mmd.clams2.0*)
