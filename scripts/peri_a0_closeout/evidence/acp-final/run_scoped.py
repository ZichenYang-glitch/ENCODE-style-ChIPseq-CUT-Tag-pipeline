"""Observe this task's compile scope from outside it, including OOM termination."""
import argparse
import json
from pathlib import Path
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('--attempt', choices=['fat', 'nonlto'], required=True)
args = parser.parse_args()
root = Path('/tmp/hw-peri-acp-final-20261109')
directory = root / args.attempt
directory.mkdir(exist_ok=False)
unit = 'hw-peri-acp-final-' + args.attempt + '.scope'
command = ['systemd-run', '--user', '--scope', '--unit=' + unit,
           '--property=MemoryMax=8G', '--property=MemorySwapMax=0',
           '--property=CPUQuota=400%', '--property=TasksMax=512',
           '/usr/bin/python3', '/tmp/helixweave-peri-design-20261108/scripts/peri_a0_closeout/verify_patch.py',
           '--source', '/tmp/hw-peri-a0-close-20261109/peri',
           '--manifest', '/tmp/hw-peri-a0-close-20261109/evidence/build-testfix.json',
           '--output', str(directory / 'rust'), '--crate', 'peri-acp',
           '--pattern', 'trusted_prediction_policy_is_strict_and_default_compatible',
           '--build-jobs', '1', '--test-opt-level-zero', '--verbose']
if args.attempt == 'nonlto':
    command += ['--test-no-lto']
record = {'command': command, 'started_at_unix': time.time(), 'sample_interval_seconds': 0.2,
          'peak_observed_memory_current': 0, 'peak_observed_kernel_memory_peak': 0,
          'observed_pids': [], 'samples': 0}
(directory / 'command.json').write_text(json.dumps(record, indent=2) + '\n')
pids = set()
with (directory / 'scope.log').open('xb') as log:
    process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
    group = None
    while process.poll() is None:
        if group is None:
            query = subprocess.run(['systemctl', '--user', 'show', unit, '--property=ControlGroup', '--value'], capture_output=True, text=True)
            if query.stdout.strip():
                group = Path('/sys/fs/cgroup') / query.stdout.strip().lstrip('/')
                record['cgroup'] = str(group)
        if group is not None:
            try:
                current = int((group / 'memory.current').read_text())
                peak = int((group / 'memory.peak').read_text())
                record['peak_observed_memory_current'] = max(current, record['peak_observed_memory_current'])
                record['peak_observed_kernel_memory_peak'] = max(peak, record['peak_observed_kernel_memory_peak'])
                record['last_memory_events'] = (group / 'memory.events').read_text()
                pids.update(map(int, (group / 'cgroup.procs').read_text().split()))
                record['samples'] += 1
            except OSError:
                pass
        time.sleep(0.2)
record.update(exit_code=process.returncode, finished_at_unix=time.time(), observed_pids=sorted(pids))
record['systemd_after'] = subprocess.run(['systemctl', '--user', 'show', unit, '--property=ActiveState', '--property=SubState', '--property=Result', '--property=ControlGroup'], capture_output=True, text=True).stdout
record['observed_pids_still_present'] = [pid for pid in pids if Path('/proc', str(pid)).exists()]
record['measurement_boundary'] = 'Scope aggregate cgroup memory (includes driver, cargo, compiler/linker and charged file cache). External 200ms observation can miss final peak; successful driver report contains kernel memory.peak before scope exit. Not native runtime sizing.'
(directory / 'resources.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'attempt': args.attempt, 'exit_code': process.returncode, 'observed_kernel_peak': record['peak_observed_kernel_memory_peak'], 'systemd_after': record['systemd_after']}))
raise SystemExit(process.returncode)
