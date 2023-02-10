# Infinity Laplacian Eigenfunctions

This code produces the examples for the paper "The Infinity Laplacian eigenvalue problem: reformulation and a numerical scheme
": https://arxiv.org/abs/2004.08127

```
@misc{bozorgnia2020infinity,
    title={The infinity Laplacian eigenvalue problem: reformulation and a numerical scheme},
    author={Farid Bozorgnia and Leon Bungert and Daniel Tenbrinck},
    year={2020},
    eprint={2004.08127},
    archivePrefix={arXiv},
    primaryClass={math.NA}
}
```

It computes eigenfunctions of the infinity Laplacian, i.e., functions $u\in W^{1,\infty}_0(\Omega)$ which are viscosity solutions of the PDE

$$
\begin{align*}
\begin{cases}
\min(|\nabla u| - \lambda u, -\Delta_\infty u) = 0,\quad &u>0 ,\\
-\Delta_\infty u = 0,\quad &u=0 ,\\
\max(-|\nabla u| - \lambda u, -\Delta_\infty u) = 0,\quad &u<0.
\end{cases}
\end{align*}
$$

The file ```RUN_ME_eigenfunction.m``` computes ground states, i.e., positive solutions, on different domains. 

The file ```RUN_ME_second_eigenfunction.m``` computes second eigenfunctions on several symmetric domains. For the square, even higher eigenfunctions can be computed.

The file ```RUN_ME_comparison.m``` gives a numerical confirmation of the recently proved statement (https://arxiv.org/abs/2210.03447) that the infinity harmonic potential on a square is no ground state.
