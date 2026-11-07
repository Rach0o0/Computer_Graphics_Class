# Computer Graphics in C++ (École Polytechnique)

Four graphics assignments written from scratch in C++ during the Computer Graphics course at École Polytechnique: a path tracer, color transfer with optimal transport, a **2D free surface fluid simulation based on semi discrete optimal transport**, and mesh parameterization.

Author: **Rachid Tazi** · [Portfolio](https://rach0o0.github.io) · [LinkedIn](https://www.linkedin.com/in/rachid-tazi/)

## Overview

| Folder | Topic | Highlights |
|---|---|---|
| [`TD1`](TD1) | **Path tracer** | Recursive path tracing, Monte Carlo antialiasing, depth of field, textured `.obj` meshes, **BVH acceleration** |
| [`TD2`](TD2) | **Color transfer** | Sliced optimal transport to match the color palette of one image to another |
| [`TD3`](TD3) | **Fluid simulation** | Voronoi and power diagrams, semi discrete optimal transport solved with **L-BFGS**, free surface incompressible fluid |
| [`TD4`](TD4) | **Mesh parameterization** | Tutte embedding of a 3D mesh onto a 2D disk |

## TD1 · Path tracer

A ray tracer built step by step: ray sphere intersection and shadows, then recursive path tracing with mirror reflections and indirect lighting, Monte Carlo antialiasing and depth of field, triangle meshes loaded from `.obj` files, and finally a Bounding Volume Hierarchy.

**What the BVH changed** (same 512×512 render, 32 paths per pixel, max depth 5):

| Strategy | Render time |
|---|---|
| Naive (every triangle, every ray) | 12 min 43 s |
| Single bounding box | 2 min 37 s |
| **BVH** | **5.9 s** (about 130x faster) |

The full report is in [`TD1/Project 1 report - Rachid Tazi.pdf`](TD1).

## TD2 · Color transfer with sliced optimal transport

The pixels of two images are treated as point clouds in RGB space. At each iteration, both clouds are projected on a random direction, sorted, matched one to one, and the input pixels are moved toward their match. After enough iterations, the input image takes the color palette of the model image.

## TD3 · Fluid simulation with semi discrete optimal transport

A 2D incompressible fluid with a free surface, built in three steps:

1. **Voronoi diagram**: each cell is built by clipping the unit square against bisectors (Sutherland Hodgman).
2. **Power diagram and optimal transport**: weights turn the Voronoi diagram into a power diagram, and **L-BFGS** adjusts them so that each cell reaches a target area.
3. **Free surface fluid**: an extra "air" weight leaves room for a free surface. At each time step, the Laguerre cells are recomputed, particles are pulled toward their cell centroid by a spring force, and gravity is applied.

Result: target fluid volume 0.4, measured **0.400184**. The fluid settles under gravity while each particle keeps the same area, which is exactly the incompressibility the method enforces.

A video of the simulation is available in [`TD3/output.mp4`](TD3/output.mp4). The full report is in [`TD3/REPORT ASSIGNMENT 2 - Rachid Tazi.pdf`](TD3).

## TD4 · Tutte embedding

Parameterization of a 3D mesh (`goethe.obj`): the boundary vertices are fixed on a circle, and each inner vertex is placed at the average of its neighbors, which gives a 2D embedding of the mesh.

## Build and run

Each folder has its own `Makefile`.

```bash
cd TD1
make
./render
```

TD3 needs the [libLBFGS](https://github.com/chokkan/liblbfgs) library installed in `/usr/local/lib`:

```bash
cd TD3
make
./geometry
```
