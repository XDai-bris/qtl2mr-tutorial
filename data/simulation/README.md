# How the workshop data were simulated

This folder holds the code that made the three data files of the workshop:
`data/cohort1.rds`, `data/cohortA.rds` and `data/cohortB.rds`. You do not need
it for the workshop. It is here so that anyone can see exactly how the
simulated mice and their traits were made, and why the right answers are what
they are.

| File | What it is |
|---|---|
| `simulate_world.R` | The main script: settings, then one numbered step per stage |
| `simulation_functions.R` | The helper functions, each with a short explanation |
| `truth/true_qtl.csv` | Where the 14 QTL are and what kind of founder pattern each has |
| `truth/true_founder_effects.csv` | The true effect of each founder at each QTL |
| `truth/attempt_log.csv` | The design checks of every simulated world that was tried |

## The truth built into the simulation

```
                 U  (hidden; never measured)
               /   \
        a     /     \     b
             v       v
  10 X-QTL → X ─────→ Y ← 4 Y-QTL
               θ = 0.5
```

| Quantity | Value | Meaning |
|---|---|---|
| θ, X → Y | 0.5 | One unit more X causes 0.5 units more Y |
| Y → X | 0 | There is no arrow from Y to X |
| a, U → X | √0.10 ≈ 0.316 | U explains 10% of the variance of X |
| b, U → Y | 0.5 / √0.10 ≈ 1.581 | Chosen so that a × b = 0.5 |
| Var(U) | 1 | U is pure noise, unrelated to the genome |
| X-QTL | 10 × 0.08 | Each X-QTL explains 0.08 of the variance of X |
| Y-QTL | 4 × 0.45 | Each Y-QTL adds 0.45 to the variance of Y; they do not touch X |
| Polygenic background | X: 0.05; U, Y: 0 | Small genome-wide effect on X, shared through kinship (about 0.04 of variance across these mice; see step 5 of the script) |
| Noise | U: 1; X: 0.05; Y: 0.25 | Independent for every mouse |

Each QTL sits on its own chromosome. Five X-QTL split the eight founders into
two groups, as a SNP would; the other five have three or more levels. The
patterns are listed in step 2 of `simulate_world.R`.

## Who does what

**simcross** makes the mice (step 4).
- `sim_do_pedigree()` writes the breeding plan of a DO population: the eight
  founder strains are crossed into pre-CC lines, and the DO is then bred for 20
  generations of random mating between non-siblings, with 144 breeding pairs.
  The last generation has 12 pups per pair, 1,728 mice.
- `sim_from_pedigree()` simulates meiosis along that pedigree on the 19
  autosomes (real mouse chromosome lengths). Every inherited chromosome is a
  recombined copy of a parent's, so after 20 generations each chromosome is a
  mosaic of short pieces from the founders, typically all eight (a few mice
  carry only seven on a given chromosome).
- At each of 1,500 evenly spaced markers we read the two founder pieces a mouse
  carries (one from each parent). This gives exact founder probabilities: 1 for
  a founder carried twice, ½ + ½ for two different founders. In real data qtl2
  estimates these from genotypes; here we know them.

**qtl2founderSim** makes the traits (steps 5 to 8).
- `make_empirical_design()` is where simcross hands over to qtl2founderSim. It
  receives the founder probabilities of all 1,728 last-generation mice (step 4)
  and the qtl2 kinship computed from them, draws the 1,600 mice of the study
  population, and keeps their kinship for the polygenic background.
- Each QTL's founder pattern is scaled so that it explains its target variance
  in these 1,600 mice (`make_founder_effects()` in the helper file).
- `causal_spec()` writes the model: three traits U, X and Y; the arrows
  U → X, U → Y and X → Y in the matrix `B` (`B[to, from]`); U marked as latent
  (hidden); the founder effects; the polygenic and noise variances; and the
  measurement plan.
- `two_sample()` is the measurement plan. It shuffles the 1,600 mice once and
  deals them out: the first 400 have X and Y measured (Cohort 1), the next 500
  only X (Cohort A), the last 700 only Y (Cohort B).
- `simulate_multitrait()` draws U, X and Y for all 1,600 mice from this one
  model, then keeps each trait only in the mice where it is measured.

**qtl2** computes each cohort's kinship (leave one chromosome out, and overall)
and runs the genome scans (`scan1`) that come with the data, the same scans a
qtl2 user would run. **qtl2mr** is used only to pick instruments for the
design checks below.

## Why the association is about 1 while the causal effect is 0.5

Write the two traits as sums of their parts:

```
X = a·U + (10 X-QTL) + polygenic + noise
Y = 0.5·X + b·U + (4 Y-QTL) + noise
```

The variance of X:

```
Var(X) = a²·Var(U) + 10 × 0.08 + polygenic + noise
       ≈ 0.10      + 0.80     + 0.05      + 0.05
       ≈ 1.00
```

The covariance of X and Y has two parts, the causal path and the shared
cause U:

```
Cov(X, Y) = 0.5·Var(X)  +  b·Cov(U, X)
          ≈ 0.5 × 1.00  +  b·a·Var(U)
          = 0.50        +  0.50
          ≈ 1.00
```

So the regression slope of Y on X is Cov(X, Y) / Var(X) ≈ 1.00 / 1.00 = 1.0:
half of it is the effect of X, the other half is U. More mice only make this
slope more precise, around the wrong number.

For completeness, the variance of Y is

```
Var(Y) = 0.5²·Var(X) + b²·Var(U) + 2·0.5·b·Cov(X, U) + 4 × 0.45 + 0.25
       ≈ 0.25        + 2.50      + 0.50              + 1.80     + 0.25
       ≈ 5.30
```

so the correlation of X and Y is about 1.00 / √5.30 ≈ 0.43.

The simulated mice come out close to these values. The variance of X is 0.97
in the 900 mice with X measured, and the variance of Y is 4.95 in the 1,100
with Y measured. In Cohort 1 (the 400 mice with both), Var(X) = 0.91 and
Cov(X, Y) = 0.79, so the slope is 0.87 and the correlation 0.38: one sample of
400 mice scatters around the population values. (Sample variances, as R's
`var()` computes them.)

**Why MR finds 0.5.** U is simulated independently of the genome: it has no
QTL and no polygenic part. A founder allele at an X-QTL therefore changes X
without being linked to U, and it reaches Y only through X. At every X-QTL the
founder effects on Y are 0.5 times the founder effects on X, and that ratio is
what MR estimates.

## Why the effect of Y on X is 0

There is simply no arrow from Y to X: `B["X", "Y"]` is 0. qtl2founderSim
turns the arrows into total effects (direct plus indirect), and the script
stops unless they are exactly the truth. For the data in this repository:

```
          U   X Y
U 0.0000000 0.0 0
X 0.3162278 0.0 0
Y 1.7392527 0.5 0
```

Read each row as "the total effect on this trait of one unit of the column
trait". X → Y is 0.5 and Y → X is 0. U's total effect on Y,
b + 0.5·a ≈ 1.739, includes its path through X. The 4 Y-QTL change Y but not
X, which is why they are the instruments for testing Y → X.

## Three cohorts from one population

The 1,600 study mice were drawn at random from the 1,728 mice of the last DO
generation, and all three cohorts are random, non-overlapping samples from
these 1,600, made with the same model. Cohort 1 has 400 mice, Cohort A 500 and
Cohort B 700. No mouse is in two cohorts, but relatives (for example
littermates) can fall into different cohorts, as in real studies of one
colony. The sizes differ on purpose: the cohorts are separate studies, not
copies of each other. Cohort B, where only Y is measured, is the largest
because its size sets how precise two-sample MR can be; Cohort A is large
enough for strong X instruments.

## Seeds and design checks

Each attempt starts from one seed (202611025, then 202611026, and so on), which
fixes every random draw: the genotypes, the choice of the 1,600 mice, the
cohorts and the traits. Each simulated world must pass four checks before it
is used. They look only at the design, never at an MR result:

- **C1, founders.** Each founder makes up 9–16% of the genome, and has at least
  20 copies at every QTL in every cohort.
- **C2, X-QTL strength.** Under the workshop's instrument rule (LOD ≥ 3, and an
  F statistic of at least 4 for the joint test of the eight founder effects),
  at least 8 of the 10 X-QTL are found in Cohort 1 and in
  Cohort A (a selected instrument within 10 cM of the true QTL).
- **C3, the bias to show.** The plain regression slope of Y on X in Cohort 1 is
  between 0.8 and 1.2.
- **C4, Y-QTL.** At least 3 of the 4 Y-QTL are found as instruments for Y that
  are not X regions (X LOD < 3) in Cohort 1.

A world that fails is thrown away and the next seed is tried; the settings
do not change between attempts. No MR estimate is computed anywhere in this
folder.

The first edition of this workshop used cohorts of 400 / 400 / 400 mice. A
draft revision used 400 / 600 / 300; its two-sample MR estimate was 0.83 (95%
interval from just above 0.50 to 1.17), which just missed the true value of 0.5
and sat nearer the association than the truth. We then enlarged Cohort B to 700
mice for precision and set Cohort A to 500, and fixed all settings, including
the seeds above, before the draw used here.

For the data in this repository the first seed, 202611025, failed C1 (a
founder had as few as 16 copies at a QTL in a cohort), and the second seed,
202611026, passed all four checks. `truth/attempt_log.csv` lists every
attempt.

## How to re-run

You need R, a C compiler for the packages below, and these packages:

| Package | Version | Where from |
|---|---|---|
| simcross | 0.10 | CRAN: `remotes::install_version("simcross", version = "0.10")` |
| qtl2 | 0.38 (used here; the script does not check this version) | CRAN: `remotes::install_version("qtl2", version = "0.38")` |
| qtl2mr | 0.6.0 | this repository: `install.packages("package/qtl2mr_0.6.0.tar.gz", repos = NULL, type = "source")` |
| qtl2founderSim | 0.0.0.9009 | GitHub: `remotes::install_github("XDai-bris/qtl2founderSim", ref = "8dfaa5a521520442195622c1ad3120e89fac62fd")` |

Then, from the project folder:

```
Rscript data/simulation/simulate_world.R
```

It takes about 2 minutes per attempt (about 4 minutes for the two attempts
behind these data) and up to about 3 GB of memory. It overwrites the three
files in `data/` and the tables in `truth/`. The script stops unless simcross 0.10, qtl2founderSim 0.0.0.9009
and qtl2mr 0.6.0 are loaded, because another version may draw different
random numbers.

The data in this repository were made with R 4.4.3 on macOS. Other versions of
R or of the packages, or another computer, may give numbers that differ in the
last digits.

## What is in `truth/`

- `true_qtl.csv`: one row per QTL: `id` (QX1–QX10 for X, QY1–QY4 for Y),
  `trait`, `chr`, `frac` (relative position on the chromosome), `type` (the
  kind of founder pattern), `marker` and `pos` (cM).
- `true_founder_effects.csv`: one row per QTL and founder: the founder's true
  effect on the QTL's trait, centred so the eight founders average 0.
- `attempt_log.csv`: one row per attempted world, with every design-check
  value and whether it passed.

## Simplifications

Real data are messier. Here the founder probabilities are exact (no genotyping
error), only the autosomes are used, there are no covariates such as sex or
batch, every QTL of a trait has the same size, and there is a single hidden
confounder. These choices keep the workshop focused; they are not claims about
real DO studies. One side effect of the design checks: accepting only worlds
whose X instruments look strong in Cohort 1 can nudge the one-sample MR
estimate slightly toward the association.
