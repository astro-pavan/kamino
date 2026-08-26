"""Concatenate the per-shard outputs of make_crust_compositions.jl into crust_compositions.csv.

`make_crust_compositions.jl --slice n/k` writes crust_compositions_slice{n}_of_{k}.csv, each
holding a disjoint subset of the Mg/Si axis. This merges them, sorts to the canonical
(mg_si, delta_iw) order, and REFUSES to write a table that is not a complete rectangle --
crust_composition.load_crust_table would reject it anyway, and a hole in the grid that reached
the interpolator would be silently smeared over rather than reported.

    python merge_crust_slices.py [--pattern GLOB] [--out PATH]
"""

import argparse
import glob
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def merge(pattern: str, out: str) -> int:
    paths = sorted(glob.glob(pattern))
    if not paths:
        sys.exit(f'no shard files matching {pattern}')

    frames, headers = [], []
    for p in paths:
        with open(p) as fh:
            headers.append([ln for ln in fh if ln.startswith('#')])
        frames.append(pd.read_csv(p, comment='#'))

    # Every shard must have been generated with the same closure and melting path, or the table
    # is a chimera. The header comments carry both.
    if len({tuple(h) for h in headers}) != 1:
        sys.exit('shards disagree on closure/path -- regenerate them all with the same settings')

    df = pd.concat(frames, ignore_index=True).sort_values(['mg_si', 'delta_iw'])
    n_mg, n_iw = df['mg_si'].nunique(), df['delta_iw'].nunique()
    if len(df) != n_mg * n_iw:
        missing = {(m, d) for m in df['mg_si'].unique() for d in df['delta_iw'].unique()} - \
                  set(map(tuple, df[['mg_si', 'delta_iw']].to_numpy()))
        sys.exit(f'incomplete grid: {len(df)} rows for {n_mg} x {n_iw}. '
                 f'Missing {sorted(missing)[:8]}{"..." if len(missing) > 8 else ""}. '
                 f'Rerun those shards -- do not merge a table with holes.')

    with open(out, 'w') as fh:
        fh.writelines(headers[0])
        fh.write(f'# Merged from {len(paths)} shards by merge_crust_slices.py\n')
        df.to_csv(fh, index=False)

    warned = df['warnings'].notna().sum() if 'warnings' in df else 0
    print(f'Written {out}: {len(df)} rows ({n_mg} Mg/Si x {n_iw} dIW), {warned} carrying warnings.')
    if warned:
        for _, r in df[df['warnings'].notna()].iterrows():
            print(f"  Mg/Si={r['mg_si']:.2f} dIW={r['delta_iw']:+.1f}: {r['warnings']}")
    return len(df)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--pattern', default=os.path.join(HERE, 'crust_compositions_slice*.csv'))
    ap.add_argument('--out', default=os.path.join(HERE, 'crust_compositions.csv'))
    args = ap.parse_args()
    merge(args.pattern, args.out)
