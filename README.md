# 🧭 Vincenty Geodesic Toolkit (WGS-84)

A precise and lightweight Python toolkit implementing **Vincenty’s direct, inverse, and interpolation solutions** for geodesics on the **WGS-84 ellipsoid**.  
It provides sub-millimeter accuracy for computing **distance**, **bearing**, and **intermediate coordinates** between geographic points.

---

## 🧩 Introduction

In 1975, **Thaddeus Vincenty** formulated highly accurate mathematical solutions to calculate the **geodesic distance and bearings** between two latitude–longitude points on an ellipsoidal model of the Earth.  
Unlike spherical methods such as the **Haversine formula**, which assume a perfectly round Earth and yield results within about **0.3% accuracy**, Vincenty’s approach accounts for the Earth’s ellipsoidal shape—achieving precision up to **0.5 mm in distance** and **0.000015″ in bearing** on the selected ellipsoid.

The **Vincenty method** remains the standard for many geodetic applications in surveying, mapping, and navigation.  
However, the *inverse* form of Vincenty’s equations may fail to converge for nearly antipodal points (points on opposite sides of the globe), a known limitation also documented in **GeographicLib** test datasets.

---

## ⚙️ Installation

1️⃣ **Clone this repository**
```bash
git clone https://github.com/<yourusername>/vincenty-toolkit.git
cd vincenty-toolkit

---

1️⃣ **Check available options**
```python
python vincenty.py -h


⚖️ License

This project is distributed under the MIT License.
Original mathematical formulation © 1975 Thaddeus Vincenty.
Python adaptation and CLI integration by Md Rokonuzzaman, rokon.mist17@gmail.com