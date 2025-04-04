# GANGSTA+ Implementation Comparison

## Alignment Metrics

| Metric | Original | Modular | Difference |
|--------|----------|---------|------------|
| Sse Count1 | None | None | N/A |
| Sse Count2 | None | None | N/A |
| Aligned Sses | None | None | N/A |
| Aligned Residues | 93 | 93 | +0 |
| Rmsd | 16.806467 | 16.806467 | +0 |
| Score | 0.030879 | 0.030879 | +0.0000 |

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

Section not found in original implementation
Section not found in modular implementation

## Residue Alignment Section

Section not found in original implementation
Section not found in modular implementation

## Rmsd Calculation Section

```diff
  Calculating RMSD and transformation for 93 aligned pairs
  Extracted 93 valid point pairs for RMSD calculation
  Point pair 0: (1.391, -12.809, 133.038) -> (62.113, -3.246, 35.209)
  Point pair 1: (0.247, -14.637, 129.932) -> (64.963, -0.857, 34.344)
  Point pair 2: (2.147, -17.905, 129.436) -> (65.093, 1.124, 31.092)
  Point pair 3: (1.403, -21.642, 129.766) -> (67.633, 3.737, 32.195)
  Point pair 4: (0.686, -23.883, 132.749) -> (69.967, 4.855, 29.451)
  RMSD calculation result: 16.806467
  Transformation details:
    Translation: (110.626418, 66.585938, 149.047319)
    Rotation matrix:
      [0.885169, 0.304925, -0.351420]
      [0.092217, -0.855300, -0.509860]
      [-0.456039, 0.418905, -0.785205]
  Calculating alignment score:
    RMSD: 16.806467
    Aligned residues: 93
    Total residues: 93
    Coverage: 1.000000
    RMSD factor: 0.030879
    Final score: 0.030879
```

## Results Section

```diff
  Final alignment result:
    Aligned pairs: 93
    RMSD: 16.806467
    Score: 0.030879
```
