# WRF Postprocessor

## Overview

Fortran tools for extracting and calculating surface fields from WRF (Weather Research and Forecasting) NetCDF output, then writing ASCII grid and field files.

## Set up a new computer

### Required software

- **GNU Fortran (`gfortran`)**: the compiler used by `compile.sh`.
- **NetCDF-C and NetCDF-Fortran development libraries**: including `netcdf.mod`, `libnetcdf`, `libnetcdff`, and the `nf-config` utility. Installing only the NetCDF-C library is insufficient.
- **Bash** to run the build and batch scripts, and **Git** to clone this repository.

Install the compiler and NetCDF libraries from the same package manager so their compiler versions and CPU architectures are compatible. Package managers install underlying dependencies such as HDF5 automatically.

### Ubuntu / Debian

Run in a terminal:

```bash
sudo apt update
sudo apt install git build-essential gfortran libnetcdf-dev libnetcdff-dev netcdf-bin
```

`netcdf-bin` provides optional inspection tools such as `ncdump`. The [NetCDF-Fortran development package](https://packages.ubuntu.com/jammy/libnetcdff-dev) supplies the Fortran interface.

### macOS

Install Apple's Command Line Tools if they are not already installed:

```bash
xcode-select --install
```

Install [Homebrew](https://brew.sh/) and follow its printed instructions to add `brew` to your shell's `PATH`. Then run:

```bash
brew install git gcc netcdf netcdf-fortran
```

Homebrew's [GCC package](https://formulae.brew.sh/formula/gcc) includes `gfortran`; [netcdf-fortran](https://formulae.brew.sh/formula/netcdf-fortran) provides the Fortran library and `nf-config`. Apple's Clang compiler alone cannot compile this project.

### Windows

Use Ubuntu under [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/windows/wsl/install), then follow the Ubuntu instructions inside the Ubuntu terminal. Run the remaining commands there as well; the supplied scripts require Bash.

### Verify Fortran and NetCDF before building

```bash
command -v gfortran
gfortran --version
command -v nf-config
nf-config --version
nf-config --fc
nf-config --fflags
nf-config --flibs
```

The first two commands must locate a compiler and print its GNU Fortran version. `nf-config` must report the installed NetCDF-Fortran version and nonempty include/link flags. Its `--fc` output identifies the compiler used for NetCDF; check that it is compatible with the `gfortran` on your `PATH`.

If either command is missing, install the packages above and ensure their `bin` directory is on `PATH` before continuing. The project build below verifies that the compiler can actually compile and link against NetCDF.

## Download and build

```bash
git clone https://github.com/Synkevych/wrf_postprocessor.git
cd wrf_postprocessor
bash compile.sh
```

If you already copied or cloned the project, enter its directory and run `bash compile.sh`. Always rebuild on the new computer; do not reuse copied `.mod`, `.o`, or executable files.

The script obtains compiler and linker flags from `nf-config`, removes old build products, and creates:

- `extract_wrf_fields`: extracts and calculates fields.
- `validate_fields`: prints minimum and maximum values for the first input time step.

## Configure and run

Run all commands from the project directory. Create input/output directories and place your WRF NetCDF file at `input/test.nc` (or change `infile` accordingly):

```bash
mkdir -p input output
```

Edit `config.nml`, for example:

```fortran
&io_nml
  infile       = "input/test.nc"
  grid_outfile = "output/grid.dat"
  pmsl_outfile = "output/pmsl_"
  ntimes1      = 28
/
```

`grid_outfile` is the grid output path; `pmsl_outfile` is a filename prefix for files such as `output/pmsl_001.dat`. Create any parent output directories beforehand. Set `ntimes1` to a positive number of time steps; processing stops at the smaller of this value and the number available in the input file.

```bash
./validate_fields
./extract_wrf_fields
```

Both programs read `config.nml`. Validation prints field ranges for inspection; extraction writes the grid and one field file per processed time step. Existing output files at those paths are replaced.

To process every regular file in an input directory:

```bash
bash process_all_inputs.sh input output 28
```

The arguments are the input directory, output directory, and time-step limit (defaults: `input`, `output`, `28`). Keep only WRF input files in that directory. The script creates a separate output subdirectory per filename stem and restores the original `config.nml` when it exits.

## Troubleshooting setup

- **`gfortran: command not found`**: install GNU Fortran and check `PATH`. Installing only GCC's C compiler or Apple's Command Line Tools is insufficient.
- **`nf-config: command not found` or missing `netcdf.mod`**: install the NetCDF-Fortran development package and check `nf-config --fflags` for the correct include directory.
- **Incompatible `netcdf.mod`, missing libraries, or architecture errors**: use a matching compiler and NetCDF installation, then rerun `bash compile.sh`. On Apple Silicon, avoid mixing Intel packages under `/usr/local` with ARM packages under `/opt/homebrew`. Use `command -v gfortran` and `nf-config --fc` to detect conflicting installations.
- **Cannot open input or output files**: check the paths in `config.nml`, run from the project directory, and create the output directories before running.
