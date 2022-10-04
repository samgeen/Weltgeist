# Weltgeist
"Wind and emission of light with time-dependent gas evolution for interstellar structures and targets"

## Introduction

Spherical hydrodynamics for feedback from massive stars

This is a Python/Numpy module based on the code Virginia Hydrodynamics 1 designed to simulate 1D spherically symmetric flows around massive stars.


## Install

1. Set up a virtual python3 environment to put the code in and activate it

2. Make sure you have and gfortran installed (e.g. sudo apt install gfortran)

3. Clone the repository to your computer and go to the Weltgeist folder where the setup.py file is

4. Get the big stellar evolution files by calling "git lfs pull" in the Weltgeist folder

5. Run python3 setup.py build to build the code.

6. Run python3 setup.py install to install it to your python environment

7. If you change your mind and want it gone, call "pip uninstall weltgeist" (you might need to do this ouside the Weltgeist folder)

If this all works, great! Otherwise get in touch with the problem and I'll take a look.

## Quickstart

To get an idea of the code, go through the worked examples in examples/ 

They're written to take you gently through how the code works. 
The comments also invite you to play with the code and learn how changing things affects the result

## Advanced Use

Here is a rough overview of the different parts of the code:

vhone - This is the engine of Weltgeist. See http://wonka.physics.ncsu.edu/pub/VH-1/ for an overview. I've forced it to be in 1D spherical coordinates, but there's no reason it can't also be in 2D, 3D or another geometry with some work.

cooling_module - This treats gas cooling - radiative energy losses from gas particles colliding or otherwise losing energy

raytracing - This is a simple ray tracer that just propagates a ray from the centre of the simulation outwards

singlestar - This reads stellar evolution tables and allows the code to simulate "realistic" stars

hydro.py - This allows Weltgeist to access the hydrodynamic variables in VH1 and automatically calculates the correct unit conversions, etc

integrator.py - This controls the timestepping of the code, as well as saving and loading

units.py - This deals with physical unit conversions between VH1 and weltgeist

sources.py - This contains code for injecting feedback sources of different kinds at the centre of the simulation

gravity.py - Solves for gravity (note: currently not yet tested properly, may contain a order ~1 unit issue to fix)

radiation.py - This runs the raytracing module to calculate the effect of radiation on the gas

cooling.py - This implements cooling_module in weltgeist

analyticsolutions.py - This contains mathematical solutions to idealised test cases that can be used to benchmark the code

ionisedtemperatures.py - This reads tables containing photoionised gas temperatures for ionised nebulae

tester.py - This contains a simple tester that can be used to implement automatic testing

Full source documentation can be found at https://samgeen.github.io/Weltgeist

This documentation is automatically generated from the source code
