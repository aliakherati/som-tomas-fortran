# som-tomas-fortran

Fortran 77 reference implementation of the coupled SOM-TOMAS atmospheric chemistry + aerosol microphysics box model.

Three physical models compiled into one executable (`box.exe`):

- **SOM** — Statistical Oxidation Model: OH-driven oxidation of a lumped precursor (GENVOC) onto a carbon/oxygen product grid. Baseline formulation from Cappa & Wilson (2012).
- **SAPRC-14** — explicit gas-phase chemistry (Carter 2010). ~306 active species, ~1800 reactions, photolysis with temperature- and zenith-dependent rates. Integrated via a custom Young-Boris predictor-corrector (`integr2.f`).
- **TOMAS** — Two-Moment Aerosol Sectional microphysics (Adams & Seinfeld 2002; Tzivion, Feingold & Levin 1987). 26 size bins × 43 chemical species, moving-sectional condensation, Brownian coagulation, vapor + particle wall loss, water / NH₃ equilibrium.

## Build

```bash
cd src

# gfortran (default, works on macOS / Linux)
make -f simple_makfile

# Intel Fortran (optional)
make            # uses Makefile, targets ifort
```

## Run

```bash
cd src
./box.exe < input > ../outputs/<run_name>.out
```

`runme.py` is a Python driver that generates canonical input files from a parameter set. The first line of the input file is the run identifier; the Fortran uses it as a prefix for every output filename.

### Output files (per run, 9 of them)

Written to `../outputs/<run_name>.*`:

| Suffix | Contents |
|---|---|
| `.dat` | main trajectory |
| `.diag` | diagnostics |
| `_gc.dat` | gas-phase concentrations (SAPRC species, ppm, `REAL*4`) |
| `_noconc.dat` | number concentrations per size bin |
| `_aemass.dat` | aerosol mass per bin per species |
| `_saprcgc.dat` | SAPRC gas-phase state |
| `_spec.dat` | species list for header alignment |
| `_vl.dat` | vapor-loss diagnostics |
| `_tau.dat`, `_kw.dat`, `_wal.dat` | wall-loss diagnostics |

## Python/JAX ports

- [som-jax](https://github.com/aliakherati/som-jax) — SOM, differentiable
- [atmos-jax-common](https://github.com/aliakherati/atmos-jax-common) — shared infrastructure (build/run wrappers, goldens loader, tolerance-aware compare primitives)
- `saprc-jax` — planned
- `tomas-jax` — planned

## License

[MIT](LICENSE).
