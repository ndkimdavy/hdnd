# HDND

HDND is a Julia research prototype for radiation hydrodynamics, with classical and mesoscopic models.

The current implementation is operational in 1D. 2D and 3D are not implemented.

---

## Install Julia

On Linux/macOS:

```bash
curl -fsSL https://install.julialang.org | sh
```

Check the installation:

```bash
julia --version
which julia
```

---

## Prepare the project

From the root of the repository:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.status()'
```

If there is an issue with the Julia environment:

```bash
julia --project=. -e 'using Pkg, Dates; Pkg.gc(; collect_delay=Day(0))'
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.status()'
```

### MPI setup

Install MPI.jl:

```bash
julia --project=. -e 'using Pkg; Pkg.add("MPI"); Pkg.instantiate(); Pkg.precompile()'
```

Install the MPI launcher wrapper:

```bash
julia --project=. -e 'using MPI; MPI.install_mpiexecjl(force=true)'
```

Add it to the current `PATH`:

```bash
export PATH="$HOME/.julia/bin:$PATH"
```

Check the installation:

```bash
which mpiexecjl
julia --project=. -e 'using MPI; MPI.versioninfo()'
```

### Active AI closure

The current `input/config1d.txt` uses `mclosure = ai`. The file
`input/closure_info.json` selects the MLP model stored in
`input/closure_mlp_ir9.onnx`. `input/closure_tf_ir9.onnx` is retained as an
alternative model. Both inference models use ONNX IR version 9.

### Synchronization with Romeo

The scripts use the SSH alias `romeo` by default. Configure that alias, or set
`REMOTE` to your Romeo `user@host`. To push the project to `$HOME/hdnd`
on Romeo:

```bash
./tool/pushsync
```

`pushsync` uses the repository root as `LOCAL` and `hdnd/` as the remote
`TARGET` by default. Both paths can be overridden with environment variables.

To retrieve a job directory into the local `out/` directory, provide its
absolute path on Romeo:

```bash
TARGET="/gpfs/scratch/<romeo-user>/<job-name>/<job-id>/" \
  ./tool/pullsync
```

`pullsync` requires `TARGET`; its local destination can be changed with `LOCAL`.

---

## Run

HDND is currently executable in 1D.

### 1D single-process run

Set:

```text
nprocx = 1
```

in `input/config1d.txt`, then run:

```bash
julia -t auto --project=. \
  ./launch/app.jl -f input/config1d.txt
```

#### Quick check

Quick execution check without output:

```bash
julia -t auto --project=. \
  ./launch/app.jl -f <(cat <<'EOF'
ndim        = 1
nx          = 128
ny          = 1
nz          = 1

dx          = 7.8125e-6
dy          = 1.0
dz          = 1.0

nprocx      = 1
nprocy      = 1
nprocz      = 1

tsim        = 1.0e-12
courant     = 0.5

mode        = meso
gamma       = 1.6666666666666667
mu          = 39.948
tolsim      = 1.0e-3

gE          = [1.47e-3, 7.94e-3, 6.11e-2, 5.0e-1, 16.0]

srcop       = losalamos
fop         = input/opAr4g.txt
ratiokR     = [1.0, 1.0, 1.0, 1.0]
ratiokP     = [1.0, 1.0, 1.0, 1.0]

nspectral   = 1024
nangular    = 1024

rho0        = 1.0
vx0         = 0.0
vy0         = 0.0
vz0         = 0.0
Th0         = 1.160452e4

Er0         = [1.0, 1.0, 1.0, 1.0]
Frx0        = [0.0, 0.0, 0.0, 0.0]
Fry0        = [0.0, 0.0, 0.0, 0.0]
Frz0        = [0.0, 0.0, 0.0, 0.0]

hrho0       = [0.5, 0.0]
hvx0        = [0.5, 0.0]
hvy0        = [0.5, 0.0]
hvz0        = [0.5, 0.0]
hTh0        = [0.5, 0.0]

hEr0        = [0.5, 0.0]
hFrx0       = [0.5, 0.0]
hFry0       = [0.5, 0.0]
hFrz0       = [0.5, 0.0]

wall        = none

gmax        = 1.0e3

mcool       = none
hpcool      = [
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0
              ]

mclosure    = randomsearch
hpclosure   = [1.0e-12, 2.5e-10, -0.99, 0.99, 0.5, 1.0e-12, 10.0, 5.0, 1000.0, 42.0]
fclosure    = none

hprmodel    = [1.0e-6, 1.0e-3]

hsolver     = hllc
rsolver     = hll
retry       = 5

outpath     = none
outfile     = none
ioformat    = none
fps         = 1.0
EOF
)
```

### 1D MPI run

For an MPI execution, `nprocx` must match the number of MPI processes.

For example:

```text
nprocx = 4
```

in `input/config1d.txt`, then run:

```bash
mpiexecjl --project=. -n 4 \
  julia -t auto --project=. \
  ./launch/app.jl -f input/config1d.txt
```

### 1D run on Romeo

Romeo requires an explicit Slurm account. On Romeo, list the accounts
associated with your user by running:

```bash
sacctmgr show user "$USER" withassoc format=User,DefaultAccount,Account,Cluster -P
```

Use your authorized account from the `Account` column in
`tool/job0-romeo.slurm` (`#SBATCH --account=<romeo-account>`).

Submit the default Slurm job:

```bash
sbatch tool/job0-romeo.slurm
```

`XINPUT`, `XWORK` and `XOUT` can be provided through environment variables.
For example:

```bash
XINPUT="input/config1d.txt" \
XWORK="$HOME/hdnd" \
XOUT="/gpfs/scratch/$USER" \
sbatch --job-name=hdnd tool/job0-romeo.slurm
```

The current Slurm resource request is:

```text
2 days
2 x64cpu nodes
8 MPI processes per node
16 Julia threads per MPI process

16 MPI processes
256 CPU cores
```

For each job, the script copies the input configuration to
`<XOUT>/<job name>/<job ID>/config1d.txt`. It sets `nprocx` to the number of
MPI processes and `outpath` to the job's `out/` directory in that copy. The
original input configuration is not modified.

The job directory has the form:

```text
<XOUT>/<job name>/<job ID>/
├── config1d.txt
└── out/
    └── ...
```

---

## Job monitoring

```bash
watch -n 1 'squeue -l -u "$USER"; echo; tail -n 10 -v tool/*.out tool/*.err'
```

---

## Numerical workflow

The coupler applies a sequential splitting:

```text
                     coupler
┌─────────────────────────────────────────────┐
│                                             │
│    ┌─────────┐      ┌───────────┐      ┌────────┐
└──▶│  hydro  │ ──▶ │ radiation │ ──▶ │ source │
     └─────────┘      └───────────┘      └────────┘
```

The current 1D execution chain is:

```text
app
 │
 ▼
Manager
 │
 ├── Domain1D / MPI
 │
 ▼
Coupler1D
 │
 ├── Hydro1D
 ├── Radiation1D
 │     └── RModel1D
 │           └── Closure1D
 └── Source1D
 │
 ▼
Monitor1D
 │
 ▼
DBMEM local
 │
 ▼
Manager.sync
 │
 ▼
DBMEM global
 │
 ▼
IIO
 │
 ▼
IOParaView1D
```

---

## Parallel execution

The current 1D implementation combines MPI domain decomposition with Julia multithreading.

`Domain1D` decomposes the global mesh along the x direction:

```text
global 1D domain
        │
        ├── MPI rank 0
        ├── MPI rank 1
        ├── ...
        └── MPI rank N-1
```

Each MPI process owns a local part of the mesh. Julia threads are used inside each process for local computations.

The MPI decomposition is controlled by:

```text
nprocx = number of MPI processes
```

MPI communication is encapsulated in `Domain1D`. The hydrodynamics and radiation solvers operate on their local domains without directly manipulating MPI communicators.

For output, `Monitor1D` stores the local fields. `Manager.sync` gathers the local data on rank 0 and reconstructs the global fields before writing them through the I/O interface.

---

## Repository structure

```text
.
├── basic
│   ├── 1d
│   │   ├── closure1d.jl
│   │   ├── coupler1d.jl
│   │   ├── domain1d.jl
│   │   ├── hydro1d.jl
│   │   ├── ioparaview1d.jl
│   │   ├── monitor1d.jl
│   │   ├── radiation1d.jl
│   │   ├── rmodel1d.jl
│   │   └── source1d.jl
│   ├── constant.jl
│   ├── linesearch.jl
│   ├── log.jl
│   ├── manager.jl
│   ├── mem.jl
│   ├── oplosalamos.jl
│   ├── randomsearch.jl
│   ├── reject.jl
│   └── riemann.jl
├── doc
├── input
│   ├── closure_info.json
│   ├── closure_mlp_ir9.onnx
│   ├── closure_tf_ir9.onnx
│   ├── config1d.txt
│   └── opAr4g.txt
├── interface
│   ├── iclosure.jl
│   ├── idomain.jl
│   ├── iio.jl
│   ├── imonitor.jl
│   ├── iopacity.jl
│   └── isolver.jl
├── launch
│   └── app.jl
├── out
│   └── hdnd
│       └── hdnd.vtkhdf
├── tool
│   ├── closure_checkpoint.pt
│   ├── closure_info.json
│   ├── closure_mlp.onnx
│   ├── closure_tf.onnx
│   ├── genclosure1d.jl
│   ├── job0-romeo.slurm
│   ├── job1.slurm
│   ├── job2.slurm
│   ├── pullsync
│   ├── pushsync
│   ├── train.html
│   └── train.ipynb
├── Manifest.toml
├── Project.toml
└── README.md
```

---

## Configuration

The main operational configuration file is `input/config1d.txt`.

| Parameter   | Format / possible values           | Role                                      |
| ----------- | ---------------------------------- | ----------------------------------------- |
| `ndim`      | `1`                                | simulation dimension                      |
| `nx`        | `128`                              | number of cells along x                   |
| `ny`        | `1`                                | not implemented                           |
| `nz`        | `1`                                | not implemented                           |
| `dx`        | `7.8125e-6`                        | cell size along x                         |
| `dy`        | `1.0`                              | not implemented                           |
| `dz`        | `1.0`                              | not implemented                           |
| `nprocx`    | positive integer                   | number of MPI processes along x           |
| `nprocy`    | `1`                                | not implemented                           |
| `nprocz`    | `1`                                | not implemented                           |
| `tsim`      | `1.0e-7`                           | final simulation time                     |
| `courant`   | `0.5`                              | CFL coefficient                           |
| `mode`      | `classic`, `meso`                  | model selector                            |
| `gamma`     | `1.6666666666666667`               | adiabatic coefficient                     |
| `mu`        | `39.948`                           | effective atomic mass                     |
| `tolsim`    | `1.0e-2`                           | global numerical tolerance                |
| `gE`        | `[E1, E2, ...]`                    | radiation group bounds                    |
| `srcop`     | `losalamos`                        | opacity source                            |
| `fop`       | path                               | opacity file                              |
| `ratiokR`   | `[r1, r2, ...]`                    | Rosseland opacity scaling factors         |
| `ratiokP`   | `[r1, r2, ...]`                    | Planck opacity scaling factors            |
| `nspectral` | positive integer                   | spectral quadrature                       |
| `nangular`  | positive integer                   | angular quadrature                        |
| `rho0`      | `1.0`                              | initial density                           |
| `vx0`       | `-7.996376e3`                      | initial velocity along x                  |
| `vy0`       | `0.0`                              | not implemented                           |
| `vz0`       | `0.0`                              | not implemented                           |
| `Th0`       | `1.160452e4`                       | initial fluid temperature                 |
| `Er0`       | one value per group                | initial radiation energy per group        |
| `Frx0`      | `[0.0, ...]`                       | initial radiation flux along x            |
| `Fry0`      | `[0.0, ...]`                       | not implemented                           |
| `Frz0`      | `[0.0, ...]`                       | not implemented                           |
| `hrho0`     | `[xcut, hrel]`                     | density heterogeneity                     |
| `hvx0`      | `[xcut, hrel]`                     | x-velocity heterogeneity                  |
| `hvy0`      | `[xcut, hrel]`                     | not implemented                           |
| `hvz0`      | `[xcut, hrel]`                     | not implemented                           |
| `hTh0`      | `[xcut, hrel]`                     | temperature heterogeneity                 |
| `hEr0`      | `[xcut, hrel]`                     | radiation energy heterogeneity            |
| `hFrx0`     | `[xcut, hrel]`                     | x-flux heterogeneity                      |
| `hFry0`     | `[xcut, hrel]`                     | not implemented                           |
| `hFrz0`     | `[xcut, hrel]`                     | not implemented                           |
| `wall`      | `none`, `left`, `right`            | boundary wall selector                    |
| `gmax`      | value greater than or equal to `1` | maximum mesoscopic source gain            |
| `mcool`     | `none`, `fplanck`, `fpower`        | thin-regime cooling selector              |
| `hpcool`    | six values per group               | power-law cooling parameters              |
| `mclosure`  | `randomsearch`, `linesearch`, `ai` | closure method selector                   |
| `hpclosure` | array                              | closure parameters                        |
| `fclosure`  | path, `none`                       | external closure model file               |
| `hprmodel`  | $[\lambda_{\mathrm{thick}} \ll dx,\ \lambda_{\mathrm{thin}} \ge nx \times dx]$ | radiation model mean-free-path thresholds |
| `hsolver`   | `rusanov`, `hll`, `hlle`, `hllc`   | hydrodynamic solver                       |
| `rsolver`   | `rusanov`, `hll`, `hlle`           | radiation solver                          |
| `retry`     | non-negative integer, default `5`  | number of retries after rejection         |
| `outpath`   | path, `none`                       | output directory                          |
| `outfile`   | name, `none`                       | output basename                           |
| `ioformat`  | `vtkhdf`, `vtr`, `none`            | output format                             |
| `fps`       | positive value                     | output frame rate                         |

To disable output:

```text
outpath = none
outfile = none
ioformat = none
```

---

## Important source and vector parameters

### `gmax`

```text
gmax = 1.0e3
```

The source coupling is explicit. In `classic` mode, the laboratory-frame source is applied directly. In `meso` mode, both the energy and momentum sources are multiplied by

$$
\min\left(\frac{c}{c_{\mathrm{eff},l}}, g_{\max}\right).
$$

The value of `gmax` must be greater than or equal to one.

### `mcool` and `hpcool`

The thin-regime cooling selector is:

```text
mcool = none
```

The available values are:

| Value      | Thin-regime source                                               |
| ---------- | ---------------------------------------------------------------- |
| `none`     | full Planck absorption and emission source                       |
| `fplanck`  | Planck emission cooling                                           |
| `fpower`   | power-law cooling                                                  |

For `fpower`, each radiation group uses six consecutive parameters:

```text
hpcool = [C0, a, b, c, d, e, ...]
```

with

$$
\Lambda_l = C_{0,l}\rho^{a_l}p^{b_l}x^{c_l}y^{d_l}z^{e_l}.
$$

Therefore, `hpcool` must contain exactly six values per radiation group. In the current argon configuration, the first three groups use the bremsstrahlung approximation `C0 = 1.3e11`, `a = 1.5`, `b = 0.5`; the fourth group is disabled.

### `hprmodel`

$$
\mathrm{hprmodel} = [\lambda_{\mathrm{thick}}, \lambda_{\mathrm{thin}}].
$$

In `meso` mode, the radiation regime is selected from the local photon mean free path $\lambda_l = 1/(\rho\kappa_{R,l})$:

* `thick` when $\lambda_l \le \lambda_{\mathrm{thick}}$;
* `m1` when $\lambda_{\mathrm{thick}} < \lambda_l < \lambda_{\mathrm{thin}}$;
* `thin` when $\lambda_l \ge \lambda_{\mathrm{thin}}$.

The thresholds are chosen such that $\lambda_{\mathrm{thick}} \ll dx$ and $\lambda_{\mathrm{thin}} \ge nx \times dx$.

In `classic` mode, the M1 closure is used independently of these thresholds.

### `hpclosure` with `mclosure = randomsearch`

```text
hpclosure = [Alphamin, Alphamax, Betamin, Betamax, shrink, epsilon, nmax, stallmax, nsamp, seed]
```

| Position | Name       | Role                                |
| -------: | ---------- | ----------------------------------- |
|        1 | `Alphamin` | lower bound for `Alpha`             |
|        2 | `Alphamax` | upper bound for `Alpha`             |
|        3 | `Betamin`  | lower bound for `Beta`              |
|        4 | `Betamax`  | upper bound for `Beta`              |
|        5 | `shrink`   | adaptive domain reduction factor    |
|        6 | `epsilon`  | convergence tolerance               |
|        7 | `nmax`     | maximum number of iterations        |
|        8 | `stallmax` | stopping criterion after stagnation |
|        9 | `nsamp`    | number of samples                   |
|       10 | `seed`     | random seed                         |

### `hpclosure` with `mclosure = linesearch`

```text
hpclosure = [Alphamin, Alphamax, Betamin, Betamax, amin, amax, epsilon, nmax, cgid, delta, nsamp, seed]
```

| Position | Name       | Role                                |
| -------: | ---------- | ----------------------------------- |
|        1 | `Alphamin` | lower bound for `Alpha`             |
|        2 | `Alphamax` | upper bound for `Alpha`             |
|        3 | `Betamin`  | lower bound for `Beta`              |
|        4 | `Betamax`  | upper bound for `Beta`              |
|        5 | `amin`     | lower bound of the Line Search step |
|        6 | `amax`     | upper bound of the Line Search step |
|        7 | `epsilon`  | convergence tolerance               |
|        8 | `nmax`     | maximum number of iterations        |
|        9 | `cgid`     | conjugate-direction choice          |
|       10 | `delta`    | finite-difference step              |
|       11 | `nsamp`    | number of starting points           |
|       12 | `seed`     | random seed                         |

---

## Environment variables

| Variable            | Possible values                     | Role                                                       |
| ------------------- | ----------------------------------- | ---------------------------------------------------------- |
| `REAL`              | `float64`, `f64`, `float32`, `f32`  | selects the floating-point type                            |
| `LOG`               | `all`, `none`, `hdnd`               | controls log output                                        |
| `REJECT`            | `all`, `none`, or filter words      | controls rejection diagnostics                             |
| `USE_LS_BRENT`      | `true`, `false`                     | selects Brent or golden-section search in `linesearch.jl`  |
| `JULIA_NUM_THREADS` | positive integer, `auto`            | sets the number of Julia threads                           |

For `tool/job0-romeo.slurm`, the following variables can also be provided externally:

| Variable | Default                 | Role                     |
| -------- | ----------------------- | ------------------------ |
| `XINPUT` | `input/config1d.txt`    | input configuration file |
| `XWORK`  | `$HOME/hdnd`            | project workspace        |
| `XOUT`   | cluster storage         | Slurm output root        |

The default `XOUT` is `/gpfs/scratch/$USER` on Romeo.

---

## Outputs

With:

```text
outpath = out/hdnd
outfile = hdnd
ioformat = vtkhdf
```

the output file is:

```text
out/hdnd/hdnd.vtkhdf
```

The `vtkhdf` format stores the simulation time series in a single file. The `vtr` format stores the time series using a `.pvd` index and the associated `.vtr` files.

Open a VTKHDF result with:

```bash
paraview out/hdnd/hdnd.vtkhdf
```

For Romeo Slurm jobs, the configuration and simulation outputs are
stored separately:

```text
<XOUT>/<job name>/<job ID>/
├── config1d.txt
└── out/
    └── ...
```

The stored `config1d.txt` is the configuration actually used by the simulation.

---

## Current state

### Currently present

* operational 1D radiation hydrodynamics
* available models: classical and mesoscopic
* available coupling: hydrodynamics / radiation / source
* MPI decomposition of the 1D domain
* hybrid MPI and Julia multithreading
* source coupling ODEs treated with an explicit scheme
* thin-regime cooling available with `fplanck` and `fpower`
* closure used in the current main configuration: `ai` with the MLP IR9 model
  selected by `input/closure_info.json`
* closure dataset and model tools: `tool/genclosure1d.jl`, `tool/train.ipynb`
* output available in `vtkhdf` and `vtr` formats
* opacities:

  * Los Alamos TOPS table for argon (`Ar4g`)
  * current opacity file: `input/opAr4g.txt`
  * source website: [Los Alamos TOPS](https://aphysics2.lanl.gov/)

### Improvements

* robust treatment to preserve the radiation constraint $|F| < c^\star E$
* robust treatment to preserve positive thermal energy $E_{\mathrm{th}} = E_h - E_k > 0$
* implicit or asymptotic treatment of stiff sources
* `linesearch`: not robust for non-Lipschitz functions

### Not implemented

* 2D simulation
* 3D simulation
* implicit source coupling
