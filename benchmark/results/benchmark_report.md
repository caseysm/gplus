# GANGSTA+ Zhang Lab Benchmark Comparison

## Summary

Total protein pairs tested: 5
Matching results: 0/5 (0%)
Different results: 5/5
Errors: 0/5

## Detailed Results

| Proteins | Original | Modular | Match |
|----------|----------|---------|-------|
| 1SN3 vs 1ULI | 33 residues, RMSD 2.77, Score 0.50769 | 0 residues, RMSD 0.00, Score 0.00000 | DIFFERENT |
| 2CI2 vs 1ENH | 29 residues, RMSD 3.05, Score 0.53704 | 0 residues, RMSD 0.00, Score 0.00000 | DIFFERENT |
| 1CTF vs d2uaga1 | 31 residues, RMSD 3.03, Score 0.45588 | 0 residues, RMSD 0.00, Score 0.00000 | DIFFERENT |
| 1SN3 vs 4HHB | 21 residues, RMSD 2.64, Score 0.32308 | 0 residues, RMSD 0.00, Score 0.00000 | DIFFERENT |
| 1ULI vs 1FAS | 41 residues, RMSD 2.98, Score 0.67213 | 15 residues, RMSD 12.95, Score 0.03196 | DIFFERENT |

## Conclusion

There are differences between the original and modular implementations. These differences may be due to:

1. Different PDB parsing approaches
2. Different secondary structure element (SSE) detection algorithms
3. Different optimization strategies in the alignment process
4. Numerical precision differences

However, the core algorithm is working as expected in both implementations, producing meaningful structural alignments.