#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (2017) Daniela Frömberg
#
# Simple test of the revreaddy software
#
from __future__ import print_function

import revreaddy as rdy
import numpy as np
from os import makedirs, path

output_path = path.join(path.dirname(path.realpath(__file__)), "..", "data")
output_path = path.join(output_path, "sim05q")
print(output_path)
if not path.exists(output_path):
    makedirs(output_path)


"""
"RevDiffusion"
"RevReactions"
"RevDiffusionRevReactions"
"FractionalDiffusion"
"""
sim = rdy.Sim("FractionalDiffusion")

sim.alpha = 0.5
sim.boxsize = 10.
sim.is_periodic = True

sim.new_type(name="S", radius=1, diffusion_constant=1.) # 0
#sim.new_type(name="P", radius=1, diffusion_constant=1.) # 1
#sim.new_type(name="C", radius=1, diffusion_constant=0.) # 2
#sim.new_type(name="E", radius=1, diffusion_constant=0.) # 3


#k_1 = 5.
#sim.new_enzymatic(name="MM", forward_type_a=0, backward_type_b=2, catalyst_type_c=3, forward_rate=k_1, backward_rate=0., reaction_distance=3.)

for i in range(5000):
	for s in [0]: # particle types: [0, 1, 2]:
		pos = np.random.random(size=3) * 10. - 5.
		sim.add_particle(init_pos=pos, particle_type_id=s)

sim.new_mean_squared_displacement(rec_period=10, filename=path.join(output_path, "msd05.dat"), particle_type_id=0)
sim.new_mean_quartic_displacement(rec_period=10, filename=path.join(output_path, "mqd05.dat"), particle_type_id=0)
sim.new_trajectory_unique(rec_period=1, filename=path.join(output_path, "traj05.h5"), clear_period=1000)

sim.run(steps=20000, timestep=0.1)

sim.write_observables_to_file()
