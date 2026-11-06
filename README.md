# SoftFrame

SoftFrame renders triangles on the CPU. A bounding volume hierarchy lets a ray skip whole groups of triangles. The renderer uses the same scene to answer two questions: how quickly can it find the closest hit, and where should it spend extra samples to reduce visible noise?

## Render the scene

Use Linux or Ubuntu in WSL on an x86 machine with Make and a C++17 compiler that supports OpenMP. Run `make test` for the geometry checks. Run `make benchmark` to render the scene and write `results/softframe.ppm`. Set OMP_NUM_THREADS to 9 to match the thread count in the saved timing.

## Read the measurements

In the saved run, nine threads traced 7.82 million coherent primary rays per second through 50,000 triangles. That timer measures intersection traversal. The separate image test compares finished renders against a reference at 1,024 samples per pixel. On this scene, adaptive allocation with Sobol sampling took 46.5% less time than fixed sampling at comparable image error. The measurements are in `results/softframe_current.json` and `results/softframe_quality.json`.

Developed with AI assistance.
