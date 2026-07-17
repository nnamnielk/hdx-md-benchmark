#!/usr/bin/env python
"""Launch all 3 proteins in parallel, each as a subprocess."""
import sys, os, subprocess

UPSIDE_HOME = os.environ.get('UPSIDE_HOME', os.path.expanduser('~/.hermes/sandbox/upside2-md'))

proteins = ['nug2b', 'ubiquitin', 'lambda_repressor']
procs = []

for name in proteins:
    env = os.environ.copy()
    env['UPSIDE_HOME'] = UPSIDE_HOME
    log = os.path.join(os.path.dirname(__file__), name, 'outputs', f'{name}_run', f'{name}.run.log')
    os.makedirs(os.path.dirname(log), exist_ok=True)
    
    p = subprocess.Popen(
        ['python', __file__, '--run-one', name],
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
    procs.append((name, p))
    print(f"Started {name} (PID {p.pid})")

print("\nWaiting for all...")
for name, p in procs:
    rc = p.wait()
    status = "✓" if rc == 0 else f"✗ (rc={rc})"
    print(f"  {status} {name}")
print("ALL DONE")
