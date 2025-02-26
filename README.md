# POKÉ

A SageMath proof-of-concept implementation of the isogeny-based public-key encryption protocol POKÉ, accompanying the research paper [POKÉ: A Compact and Efficient PKE from Higher-dimensional Isogenies](https://eprint.iacr.org/2024/624) by Andrea Basso and Luciano Maino.

The implementation relies on the SageMath implementation of (2, 2)-isogenies by Dartois, Maino, Pope, and Robert [1].


### How to run

Running the following command will produce benchmarks, as well as estimate the public key and ciphertext sizes, for all three security levels.

```
sage POKE_PKE.sage
```

<br>

---

[1]. Pierrick Dartois, Luciano Maino, Giacomo Pope, and Damien Robert. An Algorithmic Approach to (2,2)-isogenies in the Theta Model and Applications to Isogeny-based Cryptography. https://eprint.iacr.org/2023/1747
