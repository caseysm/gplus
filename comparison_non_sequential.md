# GANGSTA+ Implementation Comparison

## Alignment Metrics

| Metric | Original | Modular | Difference |
|--------|----------|---------|------------|
| Sse Count1 | 2 | 2 | +0 |
| Sse Count2 | 2 | 2 | +0 |
| Aligned Sses | 2 | 2 | +0 |
| Aligned Residues | 32 | 32 | +0 |
| Rmsd | 0.084956 | 0.084956 | +0 |
| Score | 0.029232 | 0.029232 | +0.0000 |

## Initialization Section

```diff
  Debug mode enabled
  Configuration parameters:
    Core distance: 4
    Residue distance: 8.000000
    Inversion allowed: no
    Evaluation depth: 10
    Result count: 5
    Core delta: 7
    Distance max: 11.000000
    Rescale: 5.000000
```

## Sse Alignment Section

```diff
  Step 1: Aligning secondary structure elements
  SSE alignment starting
    Molecule 1 has 2 SSEs
    Molecule 2 has 2 SSEs
    Core delta (max SSE length difference): 7
    Trying to match SSE 0 from molecule 1:
      Type: Helix, Length: 13
      SSE 0 from molecule 2: Type: Helix, Length: 13, Score: 1.000000
        New best match: SSE 0 with score 1.000000
      SSE 1 from molecule 2 has different type, skipping
    Matched SSE 0 from molecule 1 with SSE 0 from molecule 2, score: 1.000000
    Trying to match SSE 1 from molecule 1:
      Type: Strand, Length: 6
      SSE 0 from molecule 2 has different type, skipping
      SSE 1 from molecule 2: Type: Strand, Length: 6, Score: 1.000000
        New best match: SSE 1 with score 1.000000
    Matched SSE 1 from molecule 1 with SSE 1 from molecule 2, score: 1.000000
  SSE alignment complete, found 2 SSE pairs
  SSE alignment results (2 pairs):
    Pair 0: (0, 0)
    Pair 1: (1, 1)
```

## Residue Alignment Section

```diff
  Step 2: Extending SSE alignment to residue level
```

## Rmsd Calculation Section

```diff
  Calculating RMSD and transformation for 19 aligned pairs
  Extracted 19 valid point pairs for RMSD calculation
  Point pair 0: (1.391, -12.809, 133.038) -> (62.113, -3.246, 35.209)
  Point pair 1: (0.247, -14.637, 129.932) -> (64.963, -0.857, 34.344)
  Point pair 2: (2.147, -17.905, 129.436) -> (65.093, 1.124, 31.092)
  Point pair 3: (1.403, -21.642, 129.766) -> (67.633, 3.737, 32.195)
  Point pair 4: (0.686, -23.883, 132.749) -> (69.967, 4.855, 29.451)
  RMSD calculation result: 11.726266
  Transformation details:
    Translation: (163.523333, 72.939917, 82.870010)
    Rotation matrix:
      [-0.654719, 0.173913, -0.735593]
      [0.654989, 0.616260, -0.437278]
      [0.377268, -0.768100, -0.517389]
  Calculating alignment score:
    RMSD: 11.726266
    Aligned residues: 19
    Total residues: 93
    Coverage: 0.204301
    RMSD factor: 0.061431
    Final score: 0.012550
```

## Results Section

```diff
  Final alignment result:
    Aligned pairs: 19
    RMSD: 11.726266
    Score: 0.012550
  Initial residue alignment results:
    Aligned count: 19
    RMSD: 11.726266
    Score: 0.012550
  Final alignment results:
    Aligned count: 32
    RMSD: 9.845669
    Score: 0.029232
```
