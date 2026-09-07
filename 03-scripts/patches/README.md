# Patches

## `RedRibbon_qsort_r_macos.patch`

Fixes a hard crash (`segfault from C stack overflow`) that makes **every**
RedRibbon C entry point unusable on macOS — including the package's own
vignette example.

**Cause.** RedRibbon bundles the `ale` C library, whose `sort_q_indirect()`
(`src/ale-1.4/src/sort.c`) and `stats_p_adjust_fdr_bh()`
(`src/ale-1.4/src/stats/base.c`) call `qsort_r()` with the **glibc** argument
order:

```c
qsort_r(base, n, size, int (*compar)(const void*, const void*, void*), void *arg);
```

The macOS/BSD `qsort_r()` takes the opposite order, and its comparator takes the
thunk *first*:

```c
qsort_r(base, n, size, void *thunk, int (*compar)(void*, const void*, const void*));
```

So on macOS the closure pointer is passed where the comparator is expected and
libc calls a *struct on the stack* as a function → `EXC_BAD_ACCESS`. Every
RedRibbon call goes through `rrho_init()` → `sort_q_indirect()`, so everything
crashes, at any list length.

The patch adds an `ALE_QSORT_R()` shim that trampolines through the BSD
signature on Apple/BSD platforms and expands to a plain `qsort_r()` elsewhere,
so Linux behaviour is unchanged.

Upstream: <https://github.com/antpiron/RedRibbon> @ `d945366a9a829dc9372c8868f56415e9b61ca849` (v1.4-1).

### Apply

```sh
git clone https://github.com/antpiron/RedRibbon.git
cd RedRibbon && git apply ../03-scripts/patches/RedRibbon_qsort_r_macos.patch
R CMD INSTALL .
```

### Known remaining platform caveat (not fixed by this patch)

RedRibbon avoids p-value underflow by accumulating in C `long double`; the paper
quotes a floor of `3.36e-4932` instead of `2.23e-308`. That relies on x86's
80-bit extended precision. On **Apple Silicon `long double` is just `double`**
(`sizeof(long double) == 8`, `LDBL_MIN == 2.2250738585072014e-308`), so RedRibbon
saturates at `-log10 p ≈ 307.7` here exactly as base-R `phyper()` does, and
`quadrants()` reports `pvalue = 0`, `log_pvalue = Inf` for large concordant
lists. See
`06-reports/exploratory_notebooks/investigate_RRHO_map_saturation.qmd`.
