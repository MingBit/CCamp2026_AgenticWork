# External tools (git-ignored)

## MALLET 2.0 (v202108)

Topic modelling backend for pycisTopic (`run_cgs_models_mallet`). Requires Java
(tested with OpenJDK 21). If the analysis host has no Java, `conda install -n scenicplus -c conda-forge openjdk`.

```bash
cd tools
curl -fLO https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
tar xzf Mallet-202108-bin.tar.gz && rm Mallet-202108-bin.tar.gz
Mallet-202108/bin/mallet   # prints the command list
```

MALLET's Java heap is set by the `MALLET_MEMORY` environment variable (e.g. `100G`;
`pycistopic topic_modeling mallet --memory` sets it). Size it to the server: it is the main
memory consumer of topic modelling.
