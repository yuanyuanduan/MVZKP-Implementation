#!/usr/bin/env python3
"""
Generate Bristol Fashion circuit files for *arithmetic* ZK in this repo.

In OneRound_SIF_Arith (mvzkp_arith.h), wires are field elements and:
  - Gate line with AND  → field multiplication (mult_mod)
  - Gate line with XOR  → field addition (add_mod)
  - NOT / INV           → treated as wire copy (see mvzkp_arith.h)

The on-disk format matches emp::BristolFashion::from_file (circuit_file.h).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def write_header(
    f,
    num_gate: int,
    num_wire: int,
    input_group_sizes: list[int],
    output_group_sizes: list[int],
) -> None:
    """BristolFashion header: num_gate num_wire, then input groups, then output groups."""
    f.write(f"{num_gate} {num_wire}\n")
    f.write(str(len(input_group_sizes)))
    for s in input_group_sizes:
        f.write(f" {s}")
    f.write("\n")
    f.write(str(len(output_group_sizes)))
    for s in output_group_sizes:
        f.write(f" {s}")
    f.write("\n\n")


def gen_mul_chain(f, num_mul: int) -> None:
    """
    Chain of field multiplications: w2 = w0*w1, w3 = w0*w2, ... (same topology as ANDgen_bool.py).
    2 inputs (wires 0,1), 1 output wire (last in chain).
    num_wire = num_mul + 2
    """
    if num_mul < 1:
        raise ValueError("num_mul must be >= 1")
    num_gate = num_mul
    num_wire = num_mul + 2
    write_header(f, num_gate, num_wire, [1, 1], [1])

    f.write("2 1 0 1 2 AND\n")
    for i in range(1, num_mul):
        f.write(f"2 1 0 {i + 1} {i + 2} AND\n")


def gen_add_chain(f, num_add: int) -> None:
    """
    Chain of field additions: w2 = w0+w1, w3 = w0+w2, ... (XOR lines).
    """
    if num_add < 1:
        raise ValueError("num_add must be >= 1")
    num_gate = num_add
    num_wire = num_add + 2
    write_header(f, num_gate, num_wire, [1, 1], [1])

    f.write("2 1 0 1 2 XOR\n")
    for i in range(1, num_add):
        f.write(f"2 1 0 {i + 1} {i + 2} XOR\n")


def gen_add_then_mul(f) -> None:
    """Small demo: out = (w0 + w1) * w2 — 1 XOR then 1 AND."""
    num_gate = 2
    num_wire = 4  # in: 0,1,2 — tmp: 3 — out could be 3 as mul output
    write_header(f, num_gate, num_wire, [1, 1, 1], [1])
    f.write("2 1 0 1 3 XOR\n")
    f.write("2 1 2 3 3 AND\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Generate Bristol Fashion circuits for arithmetic ZK.")
    p.add_argument(
        "--mode",
        choices=("mul_chain", "add_chain", "add_then_mul"),
        default="mul_chain",
        help="mul_chain: repeated field muls (AND). add_chain: repeated adds (XOR). add_then_mul: tiny demo.",
    )
    p.add_argument("--n", type=int, default=10000, help="Chain length for mul_chain / add_chain.")
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("arith_circuit.txt"),
        help="Output .txt path.",
    )
    args = p.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as out:
        if args.mode == "mul_chain":
            gen_mul_chain(out, args.n)
        elif args.mode == "add_chain":
            gen_add_chain(out, args.n)
        else:
            gen_add_then_mul(out)

    print(f"Wrote {args.output} (mode={args.mode})")


if __name__ == "__main__":
    main()
