import glob
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

RAW = f'out/{_ID}/raw'
os.makedirs(f'out/{_ID}/figs', exist_ok=True)

# Energy minimization uses steepest descent: its x-axis is the minimization
# step count, NOT a physical time, so it cannot share a picosecond timeline
# with the dynamics. Plot it on its own step axis.
em = read_xvg(f'{RAW}/em_pot.xvg')
plt.plot(em[1:, 0], em[1:, 1], label='em')
plt.xlabel('minimization step')
plt.ylabel('potential (kJ/mol)')
plt.title('Energy Minimization')
plt.legend()
plt.tight_layout()
plt.savefig(f'out/{_ID}/figs/minimization.png')
plt.close()


def plot_dynamics(name, ylabel, title, fname):
    """Plot a quantity over the NVT -> NPT -> production timeline.

    Equilibration (NVT, NPT) is a single chain stitched onto one ps axis. Each
    production replica is independent and starts at the end of NPT, so they are
    overlaid on a common time origin.
    """
    nvt = read_xvg(f'{RAW}/nvt_{name}.xvg')
    npt = read_xvg(f'{RAW}/npt_{name}.xvg')
    npt[:, 0] += nvt[-1, 0]
    md_start = npt[-1, 0]

    plt.plot(nvt[1:, 0], nvt[1:, 1], label='nvt')
    plt.plot(npt[1:, 0], npt[1:, 1], label='npt')

    for path in sorted(glob.glob(f'{RAW}/md_lang_r*_{name}.xvg')):
        rep = os.path.basename(path).split('_')[2]  # e.g. r1
        md = read_xvg(path)
        md[:, 0] += md_start
        plt.plot(md[1:, 0], md[1:, 1], label=f'md {rep}', alpha=0.7)

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
