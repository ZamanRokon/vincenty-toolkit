#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vincenty Geodesic Toolkit
Implements Direct, Inverse, and Interpolation on the WGS-84 ellipsoid.
Author: MD Rokonuzzaman, rokon.mist17@gmail.com
"""

import math
import argparse
import pandas as pd
from pathlib import Path

# -------------------- WGS-84 CONSTANTS --------------------
a = 6378137.0
f = 1 / 298.257223563
b = a * (1 - f)
EPS = 1e-12
MAX_ITER = 200

# =========================================================
# Helper functions
# =========================================================
def to_radians(deg): return math.radians(deg)
def to_degrees(rad): return math.degrees(rad)
def normalize_lon(lon): return (lon + 180) % 360 - 180
def normalize_azimuth(az): return az % 360

# =========================================================
# Vincenty Inverse (Distance & Bearings)
# =========================================================
def vincenty_inverse(lat1, lon1, lat2, lon2):
    φ1, φ2 = map(to_radians, [lat1, lat2])
    L = to_radians(lon2 - lon1)
    U1 = math.atan((1 - f) * math.tan(φ1))
    U2 = math.atan((1 - f) * math.tan(φ2))
    sinU1, cosU1 = math.sin(U1), math.cos(U1)
    sinU2, cosU2 = math.sin(U2), math.cos(U2)
    λ = L
    for _ in range(MAX_ITER):
        sinλ = math.sin(λ)
        cosλ = math.cos(λ)
        sinσ = math.sqrt((cosU2 * sinλ)**2 +
                         (cosU1 * sinU2 - sinU1 * cosU2 * cosλ)**2)
        if sinσ == 0:
            return 0, 0, 0
        cosσ = sinU1 * sinU2 + cosU1 * cosU2 * cosλ
        σ = math.atan2(sinσ, cosσ)
        sinα = cosU1 * cosU2 * sinλ / sinσ
        cos2α = 1 - sinα**2
        cos2σm = 0 if cos2α == 0 else cosσ - 2 * sinU1 * sinU2 / cos2α
        C = f / 16 * cos2α * (4 + f * (4 - 3 * cos2α))
        λ_prev = λ
        λ = L + (1 - C) * f * sinα * (
            σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2 * cos2σm**2))
        )
        if abs(λ - λ_prev) < EPS:
            break
    else:
        raise ValueError("Vincenty inverse formula failed to converge")

    u2 = cos2α * (a**2 - b**2) / (b**2)
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    Δσ = B * sinσ * (
        cos2σm + B / 4 * (cosσ * (-1 + 2 * cos2σm**2)
        - B / 6 * cos2σm * (-3 + 4 * sinσ**2) * (-3 + 4 * cos2σm**2))
    )
    s = b * A * (σ - Δσ)
    α1 = math.atan2(cosU2 * math.sin(λ),
                    cosU1 * sinU2 - sinU1 * cosU2 * math.cos(λ))
    α2 = math.atan2(cosU1 * math.sin(λ),
                    -sinU1 * cosU2 + cosU1 * sinU2 * math.cos(λ))
    return s, normalize_azimuth(to_degrees(α1)), normalize_azimuth(to_degrees(α2))

# =========================================================
# Vincenty Direct (Destination point)
# =========================================================
def vincenty_direct(lat1, lon1, α1_deg, s):
    φ1 = to_radians(lat1)
    λ1 = to_radians(lon1)
    α1 = to_radians(α1_deg)
    sinα1, cosα1 = math.sin(α1), math.cos(α1)
    tanU1 = (1 - f) * math.tan(φ1)
    cosU1 = 1 / math.sqrt(1 + tanU1**2)
    sinU1 = tanU1 * cosU1
    σ1 = math.atan2(tanU1, cosα1)
    sinα = cosU1 * sinα1
    cos2α = 1 - sinα**2
    u2 = cos2α * (a**2 - b**2) / b**2
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    σ = s / (b * A)
    for _ in range(MAX_ITER):
        cos2σm = math.cos(2 * σ1 + σ)
        sinσ = math.sin(σ)
        cosσ = math.cos(σ)
        Δσ = B * sinσ * (
            cos2σm + B / 4 * (cosσ * (-1 + 2 * cos2σm**2)
            - B / 6 * cos2σm * (-3 + 4 * sinσ**2) * (-3 + 4 * cos2σm**2))
        )
        σ_prev = σ
        σ = s / (b * A) + Δσ
        if abs(σ - σ_prev) < EPS:
            break
    else:
        raise ValueError("Vincenty direct formula failed to converge")
    sinσ, cosσ = math.sin(σ), math.cos(σ)
    φ2 = math.atan2(
        sinU1 * cosσ + cosU1 * sinσ * cosα1,
        (1 - f) * math.sqrt(sinα**2 + (sinU1 * sinσ - cosU1 * cosσ * cosα1)**2)
    )
    λ = math.atan2(
        sinσ * sinα1,
        cosU1 * cosσ - sinU1 * sinσ * cosα1
    )
    C = f / 16 * cos2α * (4 + f * (4 - 3 * cos2α))
    L = λ - (1 - C) * f * sinα * (
        σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2 * cos2σm**2))
    )
    λ2 = λ1 + L
    α2 = math.atan2(sinα, - (sinU1 * sinσ - cosU1 * cosσ * cosα1))
    return to_degrees(φ2), normalize_lon(to_degrees(λ2)), normalize_azimuth(to_degrees(α2))

# =========================================================
# Interpolation along the geodesic
# =========================================================
def vincenty_interpolate(lat1, lon1, lat2, lon2, csv_path):
    # Ensure relative path uses /examples/
    csv_file = Path(csv_path)
    if not csv_file.is_absolute():
        csv_file = Path(__file__).resolve().parent / "examples" / csv_file
    s_total, α1, α2 = vincenty_inverse(lat1, lon1, lat2, lon2)
    df = pd.read_csv(csv_file)
    print(f"Total geodesic distance: {s_total:.3f} m, Initial bearing: {α1:.4f}°")
    results = []
    for _, row in df.iterrows():
        s_i = float(row['distance'])
        φi, λi, _ = vincenty_direct(lat1, lon1, α1, s_i)
        results.append({
            "serial no": row['serial no'],
            "distance": s_i,
            "latitude": φi,
            "longitude": λi
        })
    out_file = csv_file.with_name("interpolated_points.csv")
    pd.DataFrame(results).to_csv(out_file, index=False)
    print(f"✅ Interpolated points saved: {out_file}")

# =========================================================
# Command-line Interface
# =========================================================
def main():
    parser = argparse.ArgumentParser(
        prog="vincenty.py",
        description="Vincenty Geodesic Toolkit — compute distance, destination, and interpolation on the WGS-84 ellipsoid.",
        epilog=(
            "Examples:\n"
            "  python vincenty.py -distance -startpoint=23.776939,97.724721 -endpoint=24.374530,84.144159\n"
            "  python vincenty.py -destination -startpoint=23.776939,97.724721 -dist=1500 -bearing=45\n"
            "  python vincenty.py -interpolate -startpoint=23.776939,97.724721 -endpoint=24.374530,84.144159 -points=sample_points.csv\n"
            "\nNote: distances are in meters, bearings in degrees."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument('-distance', action='store_true', help='Compute geodesic distance and bearings between two points')
    parser.add_argument('-destination', action='store_true', help='Compute destination point from a start point, bearing, and distance')
    parser.add_argument('-interpolate', action='store_true', help='Compute intermediate points along the geodesic path')
    parser.add_argument('-startpoint', type=str, help='Start point as "lat,lon" (e.g., 23.7769,97.7247)')
    parser.add_argument('-endpoint', type=str, help='End point as "lat,lon" (for distance or interpolation modes)')
    parser.add_argument('-dist', type=float, help='Distance in meters (for destination mode)')
    parser.add_argument('-bearing', type=float, help='Initial bearing in degrees (for destination mode)')
    parser.add_argument('-points', type=str, help='CSV filename (inside /examples/) for interpolation mode')

    args = parser.parse_args()

    if args.distance:
        lat1, lon1 = map(float, args.startpoint.split(','))
        lat2, lon2 = map(float, args.endpoint.split(','))
        s, α1, α2 = vincenty_inverse(lat1, lon1, lat2, lon2)
        print(f"Distance: {s:.3f} m\nInitial bearing: {α1:.4f}°\nFinal bearing: {α2:.4f}°")

    elif args.destination:
        lat1, lon1 = map(float, args.startpoint.split(','))
        φ2, λ2, α2 = vincenty_direct(lat1, lon1, args.bearing, args.dist)
        print(f"Destination:\nLatitude: {φ2:.8f}°\nLongitude: {λ2:.8f}°\nReverse bearing: {α2:.4f}°")

    elif args.interpolate:
        lat1, lon1 = map(float, args.startpoint.split(','))
        lat2, lon2 = map(float, args.endpoint.split(','))
        vincenty_interpolate(lat1, lon1, lat2, lon2, args.points)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
