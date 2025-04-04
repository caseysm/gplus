# GANGSTA+
## Non-Sequential Protein Structure Alignment Algorithm
### Please cite: https://onlinelibrary.wiley.com/doi/abs/10.1110/ps.035469.108
### Arbeitsgruppe Knapp - Freie Universität Berlin

https://en.wikipedia.org/wiki/MIT_License

## Building GANGSTA+

```bash
# Configure and build
cmake -B build && cmake --build build

# Run example
./build/gplus d2uaga1.pdb d1gkub1.pdb

# Clean build
rm -rf build/
```