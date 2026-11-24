# L39-aero-model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24497160.svg)](https://doi.org/10.5281/zenodo.24497160)

Six-degree-of-freedom flight-dynamics model of the **L-39** jet trainer with **polynomial approximation of the aerodynamic characteristics**, an engine model and navigation equations for orthodromic flight on the WGS-84 ellipsoid, with the diurnal rotation of the Earth taken into account.

The model extends the simplified model of aircraft motion, in which the aircraft is treated as an orthogonal basis, to a full flight-dynamics model. The nonlinear aerodynamic characteristics, traditionally stored as tables in the memory of the onboard controller, are replaced by compact polynomials in the Mach number and in the lift coefficient. They need little memory, are evaluated with a fixed number of operations (Horner's method) and have continuous derivatives of all orders.

## How to cite

The model and all results are described in the article below. **If you use this code, its results or the polynomial aerodynamic model in your work, please cite:**

> Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with polynomial approximation of the aerodynamic characteristics. *Aerospace Research in Bulgaria*, XX (submitted).

The code itself (version 1.2.0) is archived in Zenodo: https://doi.org/10.5281/zenodo.24497160

```bibtex
@article{Madzharov2026L39,
  author  = {Madzharov, Anastas},
  title   = {A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
             polynomial approximation of the aerodynamic characteristics},
  journal = {Aerospace Research in Bulgaria},
  volume  = {XX},
  year    = {2026},
  note    = {Submitted}
}
```

GitHub also shows the reference under **“Cite this repository”** (file `CITATION.cff`). Every script prints the reference after it runs and with `-h` / `--help`; `python citation.py` prints it together with the BibTeX entry.

## Files

| File | Contents |
|---|---|
| `l39.py` | Aerodynamics (polynomials), engine, standard atmosphere, gravity, equations of motion `rhs(x, u)`, Runge–Kutta step `rk4` |
| `trim.py`, `common.py` | Trimming for level flight (`trim`, `trim2` – Newton's method) and simulation `sim` (Runge–Kutta 4, step 0.01 s) |
| `trim_table.py` | Trim values for level flight → `trim.json` (Table 2 of the article) |
| `step.py` | Step deflection of the elevator by −0.3° → `step.npz` (120 s) and `long.npz` (800 s with fixed throttle and with the throttle law for a new altitude of 1500 m) (Fig. 6 of the article) |
| `figs.py`, `figs_en.py` | Figures 4–6 of the article with Bulgarian / English labels |
| `validation.py` | Comparison with independent L-39 flight-manual data (Table 3 of the article) |
| `numerics.py` | Time-step refinement, eigenvalues of the linearized model, long-term step response |
| `citation.py` | The reference to be cited |
| `trim.json`, `step.npz`, `long.npz` | Results used in the article |

State vector `x = [V, α, β, ωx, ωy, ωz, ψ, γ, ϑ, n1, B, L, h]`, control vector `u = [throttle, aileron, rudder, elevator]` (angles in radians).

## Usage

Requires Python 3.9+ with NumPy, SciPy and Matplotlib:

```
pip install -r requirements.txt
python trim_table.py      # Table 2
python step.py            # step response, saves step.npz
python figs_en.py         # Figures 4–6 (English labels)
python validation.py      # Table 3: comparison with flight-manual data
python numerics.py        # numerical verification (about 2 minutes)
python step.py --help     # help and citation
```

Example in your own code:

```python
from common import trim2, sim, d2r
x0, u0 = trim2(190., 1000.)                    # V = 190 m/s, H = 1000 m
def u(t, x):
    v = u0.copy()
    if t >= 5: v[3] = u0[3] - 0.3 * d2r        # elevator step
    return v
t, X, U = sim(x0, u, 120.)
```

## Source of the data

The aerodynamic and engine characteristics follow the manual of N. Patev (2004), *Ръководство по авиационни тренажори: Моделиране на движението на самолет L-39* [Manual of flight simulators: Modelling the motion of the L-39 aircraft], NVU "V. Levski", Faculty of Aviation, Dolna Mitropolia, ISBN 954-713-066-8. According to Patev (2004, p. 14), the relations for the aerodynamic coefficients were derived from the aerodynamic characteristics of the L-39 design specification and corrected with the results of wind-tunnel and flight tests; they are part of the technical documentation of the TL-39 flight simulator. The engine model follows pp. 46–54 of the same manual. The derivatives with respect to the control-surface deflections are taken per degree.

## Comparison with flight-manual data (`validation.py`)

| Quantity | Conditions | Model | Flight manual |
|---|---|---|---|
| Static thrust, kN | H = 0 | 17.2 | 16.9 (T.O. 1T-L39C-1) |
| Maximum level-flight speed (true), km/h | H = 0 / 5000 / 10 000 m, 4000 kg | 727 / 822 / 685 | 702 / 757 / 737 (MO SSSR, 1988) |
| Maximum rate of climb, m/s | H = 0 / 6000 / 10 000 m, 4300 kg | 22.3 / 20.2 / 1.2 | 22.0 / 10.8 / 3.4 (MO SSSR, 1988) |
| Stall speed, flaps up, km/h | 3700 / 4500 kg | 182 / 200 | 172 / 190 (T.O. 1T-L39C-1) |
| Mean roll rate, full aileron, deg/s | 380 km/h | 128 | 100–150 (T.O. 1T-L39C-1) |

Near sea level the model agrees with the flight manuals within 6 %; at high altitudes the thrust model needs refinement.

## Notes and limitations

- The attitude is integrated with Euler angles, which are singular at a pitch angle of ±90°; the model is intended for flight without vertical aerobatic manoeuvres.
- Beyond the critical angle of attack the lift coefficient is switched to a prescribed value; this post-stall representation has not been validated.
- At high altitudes and low speeds (e.g. H = 5000 m, V ≤ 160 m/s) the idle thrust of the engine model exceeds the drag; the thrust approximation in the idle region needs refinement.
- Version 1.2.0 (November 2026): the scale of the thrust function F_p4(n1) was corrected to 640 kgf according to Patev (2004), p. 50, Fig. 45 (versions 1.0.0 and 1.1.0 used 1800 kgf, which overestimated the thrust by up to 1170 kgf; the maximum static thrust is now 17.2 kN). The trim values of the throttle and of the rotor speed (Table 2) and the step response (Fig. 6) change accordingly; the aerodynamic model is unchanged. Added `validation.py`, `numerics.py` and the long-term response with a throttle law in `step.py`; the absolute angular velocity is denoted ω^a in the article; new article title.
- Version 1.1.0 (November 2026): the lateral channel, equations (19) and (20) and the polynomial F_mx5 of Table 1 of the article, was corrected according to Patev (2004), pp. 29–32 (the coefficients of F_mx5 there are in units of 10⁻²). The roll damping coefficient −0.20 is adjusted to the roll performance in the flight manual, while −1.6 from the manual gives 19 °/s.

## Author

Anastas Madzharov, Unmanned Robotic Systems Laboratory, Institute of Robotics, Bulgarian Academy of Sciences, Sofia, ORCID [0000-0002-6282-617X](https://orcid.org/0000-0002-6282-617X)

---

## На български

Модел на динамиката на полета с шест степени на свобода на учебно-тренировъчния самолет L-39 с полиномиална апроксимация на аеродинамичните характеристики, модел на двигателя и навигационни уравнения за полет по ортодромия върху елипсоида WGS-84.

**При използване на програмата, на резултатите от нея или на полиномния аеродинамичен модел, моля, цитирайте:**

> Маджаров, А. (2026). Модел на динамиката на полета с шест степени на свобода на самолет L-39 с полиномна апроксимация на аеродинамичните характеристики. *Aerospace Research in Bulgaria*, XX (под печат).

## License

MIT – see `LICENSE`. The license allows free use, provided the copyright notice is kept; citation of the article above is requested.
