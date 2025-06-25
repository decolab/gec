# Generative Effective Connectivity (GEC) Framework

**Matlab implementation of the GEC framework** developed for:

> *Kringelbach, Sanz Perl, Tagliazucchi, and Deco.*\
> “Towards a naturalistic neuroscience: mechanisms underlying the flattening of brain hierarchical organisation in movie‑watching compared to rest and task.”

This repository provides MATLAB code to analyze fMRI data from the Human Connectome project collected during movie-watching and resting-state. The code implements the **GCAT (Generative Connectivity of the Arrow of Time)** framework to study how brain hierarchical organization is modulated by naturalistic stimuli.

---

## Purpose & Scope

- Uses a hierarchy-sensitive metric: Generative Connectivity of the Arrow of Time (GCAT).
- Compares hierarchical organization across movie-watching and rest paradigms.
- Supports naturalistic neuroscience.
- Designed to be compatible with HCP 3T and 7T datasets.

---

## Repository Structure

- `*.m` scripts: Core MATLAB functions and analysis scripts.
> Note: Scripts are self-contained and can be adapted to different datasets with similar structure.

---

## Requirements

- MATLAB (version ≥ 2018a recommended)
- Human Connectome Project fMRI datasets (movie-watching, resting-state)


---

## Usage Guide

1. Download the preprocessed HCP data as needed (e.g., convert to MATLAB-compatible time series).
2. Load the time series data into the workspace.
3. Run the matlab scripts.
4. Visualize and interpret the GCAT metric.
5. Compare conditions statistically.

---

## Data Source

This code is designed to work with the **Human Connectome Project (HCP)** datasets:\
[https://www.humanconnectome.org/study/hcp-young-adult](https://www.humanconnectome.org/study/hcp-young-adult)

Please cite HCP appropriately if using their data in publications.

---

## Citation

If you use this code, please cite:

Kringelbach, M. L., Sanz Perl, Y., Tagliazucchi, E., & Deco, G. (2023). Towards a naturalistic neuroscience: mechanisms underlying the flattening of brain hierarchical organisation in movie-watching compared to rest and task. *Science Advances*.

---

## License

This project is licensed under a Creative Commons Attribution-NonCommercial 4.0 International Public License.

---

## Further Information

For questions or contributions, please contact gustavo.deco@upf.edu

Full documentation and updates will be provided at:\
[https://github.com/decolab/gec/blob/main/README.md](https://github.com/decolab/gec/blob/main/README.md)

