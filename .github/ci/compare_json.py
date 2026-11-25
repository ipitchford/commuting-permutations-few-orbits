# SPDX-License-Identifier: MIT
"""Compare a fresh results JSON with the recorded one. Policy (declared before tagging; also in README.md, ENVIRONMENT.txt,
REPLAY_RECEIPT.md and CI): the key elapsed_seconds (run time) is ignored; every other value must be identical, except that
on platforms other than the release platform, decimal strings and floats produced by mpmath may differ in the last places:
they must agree to 10 significant digits. Integers, booleans, rationals written as strings 'a/b', and all text identical.
usage: compare_json.py NEW RECORDED [--exact]"""
import json, re, sys
NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
exact = '--exact' in sys.argv
new, old = (json.load(open(p, encoding='utf-8')) for p in sys.argv[1:3])
bad, stats = [], {'identical': 0, 'close': 0}
def close(x, y):
    if x == y: return True
    if exact: return False
    if re.fullmatch(r'[-+]?\d+', str(x)) or re.fullmatch(r'[-+]?\d+', str(y)): return False
    fx, fy = float(x), float(y)
    return abs(fx - fy) <= 1e-10 * max(abs(fx), abs(fy))
def walk(a, b, path):
    if isinstance(a, dict) and isinstance(b, dict):
        ka, kb = set(a) - {'elapsed_seconds'}, set(b) - {'elapsed_seconds'}
        if ka != kb: bad.append(f'{path}: keys differ {sorted(ka ^ kb)}'); return
        for k in sorted(ka): walk(a[k], b[k], f'{path}.{k}')
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): bad.append(f'{path}: length {len(a)} vs {len(b)}'); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, f'{path}[{i}]')
    elif isinstance(a, str) and isinstance(b, str):
        if a == b: stats['identical'] += 1; return
        pa, pb = NUM.split(a), NUM.split(b); na, nb = NUM.findall(a), NUM.findall(b)
        if pa != pb or len(na) != len(nb) or not all(close(x, y) for x, y in zip(na, nb)): bad.append(f'{path}: {a!r} vs {b!r}')
        else: stats['close'] += 1
    elif isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        if a is not b and a != b: bad.append(f'{path}: {a!r} vs {b!r}')
        else: stats['identical'] += 1
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if (isinstance(a, int) and isinstance(b, int)) and a != b: bad.append(f'{path}: {a} vs {b}')
        elif not close(a, b): bad.append(f'{path}: {a} vs {b}')
        else: stats['identical' if a == b else 'close'] += 1
    else:
        bad.append(f'{path}: type {type(a).__name__} vs {type(b).__name__}')
walk(new, old, '$')
if bad:
    print('DIFFERENCES in', sys.argv[2]); [print(' ', m) for m in bad[:20]]; sys.exit(1)
print(f"OK {sys.argv[2]}: {stats['identical']} values identical, {stats['close']} within 10 significant digits (elapsed_seconds ignored)")
