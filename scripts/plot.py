import os

import matplotlib.pyplot as plt
import numpy as np


def read_xvg(filename):
    """Read a GROMACS .xvg file, skipping the grace/comment header.

    The header length is not fixed (it varies with the number of legends and
    title lines), so we filter by prefix rather than slicing a fixed count.
    """
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in ("#", "@", "&"):
                continue
            rows.append(line.split())
    return np.array(rows, dtype=float)


_ID = os.environ.get("ID")
print(_ID)

os.makedirs(f'out/{_ID}/figs', exist_ok=True)

# Energy minimization uses steepest descent: its x-axis is the minimization
# step count, NOT a physical time, so it cannot share a picosecond timeline
# with the dynamics. Plot it on its own step axis.
em = read_xvg(f'out/{_ID}/raw/em_pot.xvg')
plt.plot(em[1:, 0], em[1:, 1], label='em')
plt.xlabel('minimization step')
plt.ylabel('potential (kJ/mol)')
plt.title('Energy Minimization')
plt.legend()
plt.tight_layout()
plt.savefig(f'out/{_ID}/figs/minimization.png')
plt.close()


def stitch(stages):
    """Concatenate per-stage (time, value) arrays onto one continuous
    picosecond timeline, offsetting each stage by the previous stage's end.
    The dynamics timeline starts at t=0 (energy minimization is excluded)."""
    offset = 0.0
    for arr in stages.values():
        arr[:, 0] += offset
        offset = arr[-1, 0]
    return stages


def plot_dynamics(name, ylabel, title, fname):
    stages = stitch({
        'nvt': read_xvg(f'out/{_ID}/raw/nvt_{name}.xvg'),
        'npt': read_xvg(f'out/{_ID}/raw/npt_{name}.xvg'),
        'md': read_xvg(f'out/{_ID}/raw/md_lang_{name}.xvg'),
    })
    for label, arr in stages.items():
        plt.plot(arr[1:, 0], arr[1:, 1], label=label)
    plt.xlabel('time (ps)')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'out/{_ID}/figs/{fname}')
    plt.close()


plot_dynamics('pot', 'potential (kJ/mol)', 'Potential Energy', 'potential.png')
plot_dynamics('etot', 'total energy (kJ/mol)', 'Total Energy', 'total_energy.png')
plot_dynamics('temp', 'temperature (K)', 'Temperature', 'temperature.png')
