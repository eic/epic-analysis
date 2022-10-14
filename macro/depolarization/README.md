# Depolarization cross check

- `analysis_depolarization.C`
  - expects config file `datarec/delphes.${energy}.config` for a specific
    `$energy`
- `postprocess_depolarization.C` to draw depolarization vs. Q2 profiles
- `combine_depolarization_plots.rb` will combine depolarization profiles
  for various beam energies (depends on the `pycall` gem)
- `draw_Q2vsX.rb` draws plots of depolarization vs. (x,Q2)
