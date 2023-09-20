import os

import matplotlib.pyplot as plt
import numpy as np


def read_xvg(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    lines = lines[24:]
    lines = [line.split() for line in lines]
    lines = np.array(lines, dtype=float)
    return lines


_ID = os.environ.get("ID")
print(_ID)

os.makedirs(f'out/{_ID}/figs', exist_ok=True)

pot = {
    'em': read_xvg(f'out/{_ID}/raw/em_pot.xvg'),
    'nvt': read_xvg(f'out/{_ID}/raw/nvt_pot.xvg'),
    'npt': read_xvg(f'out/{_ID}/raw/npt_pot.xvg'),
    'md': read_xvg(f'out/{_ID}/raw/md_lang_pot.xvg'),
}

pot['nvt'][:, 0] += pot['em'][-1, 0]
pot['npt'][:, 0] += pot['nvt'][-1, 0]
pot['md'][:, 0] += pot['npt'][-1, 0]

plt.plot(pot['nvt'][1:, 0], pot['nvt'][1:, 1], label='nvt')
plt.plot(pot['npt'][1:, 0], pot['npt'][1:, 1], label='npt')
plt.plot(pot['md'][1:, 0], pot['md'][1:, 1], label='md')
plt.plot(pot['em'][1:, 0], pot['em'][1:, 1], label='em')
plt.xlabel('time (ps)')
plt.ylabel('potential (kJ/mol)')
plt.title('Potential Energy')
plt.legend()
plt.tight_layout()
plt.savefig(f'out/{_ID}/figs/potential.png')
plt.close()

etot = {
    # 'em': read_xvg(f'out/{_ID}/raw/em_etot.xvg'),
    'nvt': read_xvg(f'out/{_ID}/raw/nvt_etot.xvg'),
    'npt': read_xvg(f'out/{_ID}/raw/npt_etot.xvg'),
    'md': read_xvg(f'out/{_ID}/raw/md_lang_etot.xvg'),
}

# etot['nvt'][:, 0] += etot['em'][-1, 0]
etot['npt'][:, 0] += etot['nvt'][-1, 0]
etot['md'][:, 0] += etot['npt'][-1, 0]

plt.plot(etot['nvt'][1:, 0], etot['nvt'][1:, 1], label='nvt')
plt.plot(etot['npt'][1:, 0], etot['npt'][1:, 1], label='npt')
plt.plot(etot['md'][1:, 0], etot['md'][1:, 1], label='md')
# plt.plot(etot['em'][1:, 0], etot['em'][1:, 1], label='em')
plt.xlabel('time (ps)')
plt.ylabel('total energy (kJ/mol)')
plt.title('Total Energy')
plt.legend()
plt.tight_layout()
plt.savefig(f'out/{_ID}/figs/total_energy.png')
plt.close()

temp = {
    # 'em': read_xvg(f'out/{_ID}/raw/em_temp.xvg'),
    'nvt': read_xvg(f'out/{_ID}/raw/nvt_temp.xvg'),
    'npt': read_xvg(f'out/{_ID}/raw/npt_temp.xvg'),
    'md': read_xvg(f'out/{_ID}/raw/md_lang_temp.xvg'),
}

# temp['nvt'][:, 0] += temp['em'][-1, 0]
temp['npt'][:, 0] += temp['nvt'][-1, 0]
temp['md'][:, 0] += temp['npt'][-1, 0]

plt.plot(temp['nvt'][1:, 0], temp['nvt'][1:, 1], label='nvt')
plt.plot(temp['npt'][1:, 0], temp['npt'][1:, 1], label='npt')
plt.plot(temp['md'][1:, 0], temp['md'][1:, 1], label='md')
# plt.plot(temp['em'][1:, 0], temp['em'][1:, 1], label='em')
plt.xlabel('time (ps)')
plt.ylabel('temperature (K)')
plt.title('Temperature')
plt.legend()
plt.tight_layout()
plt.savefig(f'out/{_ID}/figs/temperature.png')
plt.close()

