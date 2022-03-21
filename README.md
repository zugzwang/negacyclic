# negacyclic

Package negacyclic implements arithmetic in the negacyclic rings
`R = Z[X]/(X^N+1)`, `R_p = Z_p[X]/(X^N+1)` and `R_{pq^l} =
Z_{pq^l}[X]/(X^N+1)` for primes `p, q` such that `p ≡ q ≡ 1 mod 2N` and `N`
is a power of 2.

It defines the `negacyclic.Polynomial` object for polynomials with
arbitrary-precision coefficients, and the `negacyclic.Vector` object, for
polynomials with small (e.g. int) coefficients.

Arithmetic in `R_p` is implemented with the Number Theoretic Transform,
arithmetic in `R_{pq^l}` is implemented with the CRT and Hensel's lifting,
and arithmetic in `R` supports Karatsuba experimentally, but by default it
chooses a prime larger than the expected coefficients and uses NTT.
