# Double Shear Bolted Connection – Abaqus Modeling (P1.py)

This repository provides a fully documented Abaqus CAE Python script for simulating **double shear bolted connections**.  
The model is associated with the research titled:

> **"A novel procedure to determine yield and ultimate load and deformation capacity of Double Shear Bolted Connections"**  
> *Manuscript under preparation; not yet published.*

---

## 🔍 DOI

This code is archived and citable via Zenodo:  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16777657.svg)](https://doi.org/10.5281/zenodo.16777657)

---

## 👨‍💻 Authors

- **Md. Ibrahim Kholil**  
  Postgraduate Student, Department of Civil Engineering, Khulna University of Engineering and Technology, Khulna, Bangladesh  
  📧 engikholil@gmail.com · [ORCID: 0000-0002-3349-8496](https://orcid.org/0000-0002-3349-8496)  

- **Khondaker Sakil Ahmed**  
  Associate Professor, Department of Civil Engineering, Military Institute of Science and Technology (MIST), Dhaka, Bangladesh  
  📧 drksa@ce.mist.ac.bd · [ORCID: 0000-0001-5010-7306](https://orcid.org/0000-0001-5010-7306)  

- **Aziz Ahmed** *(Corresponding Author)*  
  Senior Lecturer, School of Civil, Mining, Environmental and Architectural Engineering, University of Wollongong, Australia  
  📧 aziza@uow.edu.au · [ORCID: 0000-0001-9707-2606](https://orcid.org/0000-0001-9707-2606)  

---

## 🚀 Features

- Fully automated Abaqus CAE modeling using Python scripting  
- Incorporates material plasticity, symmetry boundary conditions, and surface-to-surface frictional contact  
- Validated numerical modeling framework for double shear bolted joints  

---

## 🛠️ Getting Started
📖 Citation

If you use this code, please cite it as:
@software{Kholil2025_double_shear,
  author       = {Kholil, Md. Ibrahim and Ahmed, Khondaker Sakil and Ahmed, Aziz},
  title        = {Double Shear Bolted Connection – Abaqus Modeling (P1.py)},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.16777657},
  url          = {https://doi.org/10.5281/zenodo.16777657}
}


### To run the script:

```bash
abaqus cae noGUI=P1.py

