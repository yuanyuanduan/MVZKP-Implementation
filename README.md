# MVZKP Implementation

This repository contains a prototype implementation of **MVZK (Multi-Verifier Zero-Knowledge)** protocols: zero-knowledge–style computation over circuits with multiple verifiers, covering **binary** (IKNP and PCG preprocessing variants) and **arithmetic** circuits, plus multi-verifier VOLE-style building blocks in `mv-qsvole`. It builds on the [EMP-toolkit](https://github.com/emp-toolkit) (`emp-tool`, `emp-ot`, `emp-agmpc`, `emp-zk`). The top-level CMake project name is `mvzkp`.

## Repository layout

| Path | Description |
|------|-------------|
| `mvzkp/` | MVZKP protocol headers and benchmarking wrappers for binary/arithmetic settings (e.g. `mvzkp_bool_iknp.h`, `mvzkp_bool_pcg.h`, `mvzkp_arith.h`) |
| `mv-qsvole/` | Multi-verifier (extended) VOLE and related programmatic VOLE primitives |
| `circuits/` | Circuit files (Bristol Fashion, Bristol Format, etc.) |
| `test/` | Executable benchmarks and examples (link the `mvzkp` interface library) |
| `third_party/` | EMP components and OpenSSL 1.1.1w sources (built locally by CMake) |
| `script/` | Local multi-process launch scripts (`run_3`, `run_4`) |

## Environment and dependencies

On **Ubuntu 22.04 LTS** (or similar Linux), install:

```bash
sudo apt-get update
sudo apt-get install -y cmake git build-essential perl
```

The first build compiles static OpenSSL from `third_party/openssl-1.1.1w`; ensure sufficient disk space and a working toolchain.

## Build

From the repository root:

```bash
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

Artifacts appear under `build/` (e.g. `build/test/aes_iknp_test`).

## Running tests

**Run commands from the `build/` directory** so tests resolve circuit paths under `../circuits/...`.

Each test sets **`const static int nP` to `4` by default**, so use **`script/run_4`** to start four processes. If you only launch three processes with `run_3` while `nP` remains 4, runs may hang or deadlock. `run_3` and `run_4` use different base ports (`17825` and `19825`), so you can run two experiments in parallel without port clashes.

Example (binary, Bristol Fashion, AES-128, IKNP variant):

```bash
../script/run_4 ./test/aes_iknp_test
```

Other registered test binaries:

| Binary | Description |
|--------|-------------|
| `aes_iknp_test` | Binary field, `circuits/bristol_fashion/aes_128.txt`, IKNP |
| `sha_iknp_test` | Binary field, `sha256.txt`, IKNP |
| `and_iknp_test` | Binary field, `ands_bool.txt`, IKNP |
| `and_pcg_test` | Binary field, `ands_bool.txt`, PCG |
| `and_arith_test` | Arithmetic field, `arith_circuit.txt` (path resolved; run from `build/` recommended) |
| `aes_wrk_test` | Uses `bench_mpc_once` with `circuits/bristol_format/AES-non-expanded.txt` (baseline-style benchmark) |

## Changing the number of parties

1. Edit **`const static int nP`** in the relevant **`test/*.cc`** file to match the intended party count.
2. Use the matching script: **`../script/run_3`** for three parties, **`../script/run_4`** for four. For other counts, copy a script and adjust ports as needed.

## Multi-machine or large party counts

For runs across machines or beyond the default maximum party count, edit **`third_party/emp-agmpc/emp-agmpc/cmpc_config.h`** (IP list and `NetIOMP`-related limits), consistent with EMP-agmpc.

## Batch runs

The sample command in `script/batchrun.sh` may not match the current default `nP` or script choice; align `run_3` / `run_4` and output paths with your `nP` before use.

## License

See [LICENSE](LICENSE) in the repository root (MIT).
