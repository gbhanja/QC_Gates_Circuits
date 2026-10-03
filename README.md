# `Quantum Gates — Complete Walkthrough & Figure Guide`

**Environment used:** Python 3.13 (Anaconda), Qiskit **2.5.2**, qiskit-aer **0.17.2**, matplotlib 3.10.7, numpy 2.3.5.

---

## 0. The big idea of the notebook

This notebook is the **experimental counterpart**: for almost every equation in those chapter it does three things —

1. **Builds the circuit** in Qiskit and draws it with the matplotlib circuit drawer (`style={'name':'bw'}` → the black-and-white boxes you see everywhere).
2. **Visualises the state** on the Bloch sphere / state-city / Q-sphere, or plots the amplitude/phase directly with matplotlib.
3. **Measures it** — appends `measure_all()`, runs on `AerSimulator` (2048 or 4096 shots), and shows a **counts histogram**. This turns every mathematical claim into an observable experimental result (probabilities, phases made visible, truth tables, no-cloning failure, etc.).

It is deliberately organised to *section numbering*, and it ends with a mapping table and student exercises.

### How to read the rest of this document

* Each numbered section below = one notebook section, with **what the code does**, **what the printed numbers mean**, and **what each figure shows**.
* Page references are to the PDF pages (e.g. p. 11) so you can jump straight there.
* At the end there is a **figure-type index**, an **honest review of bugs/subtleties** I found by re-running the physics, and notes for re-running the notebook.

---

## 1. Front matter and setup 

**Contents**
* Title, provenance (based on `PH 643: Quantum Information and Computation course Notes`), and credits: IBM Qiskit Community tutorials, *Coding With Qiskit* ep. 4 (Gates), *Teach Me Quantum* 2018, Qiskit/qiskit-tutorials, IBM Quantum Learning.
* The `pip install qiskit qiskit-aer matplotlib numpy pylatexenc` cell, **including all its console output** (that's why it's full of "Requirement already satisfied") — typical of an exported notebook where output was never cleared.
* A numbered reference list.

**The toolkit cell ** imports `numpy`, `matplotlib`, `QuantumCircuit`, `transpile`, and from `qiskit.quantum_info`: `Statevector`, `Operator`, `DensityMatrix`, `state_fidelity`; from `qiskit.visualization`: `plot_bloch_multivector`, `plot_bloch_vector`, `plot_histogram`, `plot_state_city`, `plot_state_qsphere`; from `qiskit.circuit.library`: `QFT`, `UnitaryGate`, `MCXGate`. `AerSimulator` is wrapped in a `try/except` so the notebook still runs (minus measurement cells) if Aer is missing. It prints `Qiskit environment ready. Aer available: True`.

**Section 0 — Visualisation Toolkit**, four helper functions reused everywhere:

| Helper | What it does |
|---|---|
| `show_circuit(qc, title)` | `qc.draw('mpl', style={'name':'bw'})` + optional title → circuit picture |
| `show_bloch(sv, title)` | `plot_bloch_multivector(sv)` → one Bloch sphere per qubit |
| `measure_and_plot(qc, shots=2048, title)` | copies `qc`, calls `measure_all()`, `transpile`s, runs on Aer, returns `(counts, histogram_figure, measured_circuit)` |
| `compare_bloch(sv_in, sv_out, …)` | intended as a side-by-side before/after Bloch figure |
| `bloch_angles(sv)` (defined later, p. 15) | converts a statevector into the Bloch angles θ, φ via the density matrix: `x=2Re ρ₀₁`, `y=2Im ρ₁₀`, `z=ρ₀₀−ρ₁₁` |

> **Note:** `compare_bloch` is effectively dead code — it creates a 2-axis figure, immediately `plt.close`s it, then returns two *separate* figures. The notebook never actually uses it (the before/after comparisons are done as two `show_bloch` calls instead). It's harmless but a leftover.

---

## 2. Section 1 — The general single-qubit gate 

### 2.1 The physics being tested
Any single-qubit evolution with Hamiltonian `H = h₀·I + β n·σ` gives, up to global phase,
`U = exp(iα)·Rn(η)`, with `Rn(η) = cos(η/2)·I − i sin(η/2)(n·σ)`. `Rn(η)` **rotates the Bloch vector by η anticlockwise about axis n**.

### 2.2 What the code does
* Defines `Rn_matrix(n, eta)` — literally implements with the three Pauli matrices. This is the "theory" version.
* Defines `equal_up_to_phase(A, B)` — checks two unitaries are equal modulo a global phase (needed because Qiskit's `Rz` is the notes' `Rn` *only up to phase*).
* Builds `Rz(0.7)` in Qiskit and compares its `Operator` matrix with `Rn_matrix([0,0,1], 0.7)`.

**Printed result :** both matrices are `diag(0.939−0.343i, 0.939+0.343i)` → `Match: True`. So **Qiskit's `Rz(η)` = the notes' `Rn(η)` about n = ẑ**.

### 2.3 Z-basis vs X-basis measurement 
Circuits: `H` then `Rz(π/2)`; and the same followed by another `H` before measurement.

* Figure : **circuit diagram** `q: ─H─Rz(π/2)─` — a state-prep H (to make |+⟩) and the phase gate, plus a second circuit `─H─Rz(π/2)─H─M─` for the basis-change version.
* **Figure — "Z-basis measurement: 50/50, Rz phase invisible here" (p. 7):** two bars ≈ 1002 and 1046 out of 2048. This is the central lesson of phase gates: `Rz` multiplies a *relative phase*, and a Z-basis measurement cannot see relative phase, so the statistics are an uninformative 50/50.
* **Figure — "X-basis measurement: Rz phase now visible" (p. 8):** again ≈ 1032/1016 — still ~50/50 because π/2 happens to be a special angle; but the *point* (that the H before detection rotates phase information into Z-basis probabilities) is made by the two circuits shown on p. 9 ("Full circuit with basis-change + measurement" — H, Rz, barrier, H, measure, showing both the `meas` register and how the measurement is attached).

**Take-away:** the notebook immediately establishes the theme it returns to repeatedly — *phases are invisible unless you interfere them.*

---

## 3. Section 1.1 — Rn(η) ⇔ rotation of the Bloch vector 

Checks two of the notes' identities by watching the Bloch pointer move:

* `Rz(η)|θ,φ⟩ = |θ, φ+η⟩` 
* `Rx(η)|θ,π/2⟩ = |θ−η, π/2⟩` 

**Rz experiment (θ = π/3, φ = π/5, η = π/2):**
* Circuit : `q: ─Ry(π/3)─Rz(π/5)─` labelled *"State-preparation circuit for |θ,φ⟩"* — note it uses `Ry` to set θ and `Rz` to set φ, the standard way to reach any `|θ,φ⟩`.
* **Bloch sphere "Before Rz(eta)" (p. 12):** the whole sphere is tilted so that the north pole reads `|0⟩` (that is the state-prep `Rz(π/5)` — exactly the φ-shift the notes predict); then the notebook applies `Rz(η)`.
* Circuit : `─Ry(π/3)─Rz(π/5)─Rz(π/2)─`.
* **Bloch sphere "After Rz(eta)" (p. 15):** the arrow points at the **equator**, i.e. θ is unchanged and only the azimuth moved. This is the graphic proof.
* Printed check : `θ: 1.047 → 1.047 (unchanged)` and `φ: 0.628 → 2.199 = φ + η mod 2π`. ✔

**Rx experiment (θ = 2π/5, φ = π/2, η = π/4):**
* Circuit : `─Ry(2π/5)─Rz(π/2)─Rx(π/4)─`.
* **Bloch sphere "Before Rx(eta)" (p. 18):** arrow on the −x axis, in the y–z plane (φ = π/2) — expected from `Ry`+`Rz(π/2)`.
* **Bloch sphere "After Rx(eta)" (p. 19):** arrow now tilted up toward the north pole/​+z direction — a rotation about the x-axis, so the y-component has vanished: the state has left the equator.
* Printed check : `θ: 1.257 → 0.471 = θ − η` ✔ and `φ: 1.571 → 1.571` (unchanged, = π/2) ✔.

This is the cleanest "one row of figures = one equation" demonstration in the whole notebook.

---

## 4. Section 2 — The elementary single-qubit gates 

A single generic helper, `gate_report(name, build_fn, eq_ref, input_label)`, produces for each gate: circuit drawing → the gate **matrix** → the output **statevector** → the **Bloch sphere** of `G|input⟩`. Then, separately, a **measurement histogram**.

### 4.1 X gate 
* Circuit : `q: ─X─`; matrix `[[0,1],[1,0]]`; `X|0⟩ = |1⟩`.
* **Bloch sphere "X|0⟩ on Bloch sphere" (p. 22):** arrow points straight **down** to `|1⟩` at the south pole — the π-rotation about x.
* **Measurement circuit :** `─X─ ▧ ─M─` (gate + measurement box) and the **histogram (p. 24–25) "X|0⟩ measurement: always 1"** — one bar of height **2048/2048**. Because `X|0⟩` is a computational basis state, measurement is deterministic: the classic "single spike" signature.

### 4.2 Y gate 
* Matrix `[[0,−i],[i,0]]`, `Y|0⟩ = i|1⟩`.
* **Bloch sphere "Y|0⟩" :** again straight down to `|1⟩` — because `Y|0⟩ = i|1⟩` differs from `|1⟩` by the global phase `i`, which the Bloch sphere (correctly) cannot show. A neat, non-obvious illustration that global phase is unobservable.

### 4.3 Z gate 
* Input chosen as `|1⟩`: matrix `diag(1,−1)`, `Z|1⟩ = −|1⟩`.
* **Bloch sphere "Z|1⟩" :** looks exactly like plain `|1⟩` — the minus sign is a *global* phase for this input.
* The notebook then does the **H–Z–H trick :** `─H─Z─H─` converts the invisible phase flip into a bit flip: `HZH|0⟩ = |1⟩`.
* **Histogram "H–Z–H |0⟩: should measure |1⟩ deterministically"** — single bar at `1`, 2048/2048. This is the first *interference* demonstration in the notebook.

### 4.4 T (π/8) gate 
* Input `|1⟩`; matrix `diag(1, e^{iπ/4})`, i.e. `T|1⟩ = (0.707+0.707i)|1⟩`.
* The code independently rebuilds the notes' formula `T = e^{iπ/8}(cos(π/8)·I − i sin(π/8)·Z)` and prints **"Matches Qiskit T gate: True"** (p. 32). This is exactly the "π/8 gate" naming justification *and* the statement that `T` is an `Rz(π/4)` up to global phase.
* **Bloch sphere "T|1⟩" :** still the south pole — same "invisible phase" story as Y and Z.

### 4.5 S gate 
* `S = diag(1, i)`; check `T·T == S` → printed **True** ; circuit `─T─T─` is drawn next to `─S─` as a visual proof that the two circuits are the same operator.
* **Bloch sphere "S|1⟩" :** again the south pole (phase-only gate).

### 4.6 Hadamard gate 
* Matrix `(1/√2)[[1,1],[1,−1]]`; the code verifies the notes' Pauli form `H = (X+Z)/√2`  → **True** ).
* **Bloch sphere "H|0⟩" :** arrow points to the **+x axis** — the "half-way between x and z" π-rotation about the (x+z)/√2 axis, giving `|+⟩`.
* **Measurement circuit  and histogram  "H|0⟩ measurement: expect ~50/50 split":** bars of **2071** and **2025** out of 4096. Unlike the previous deterministic spikes, this is the genuinely random, quantum result — the notebook's first "real" quantum statistics.

---

## 5. Section 3 — Universal single-qubit gates: H and T gates 

**Physics:** repeated `H` and `T` generate a *dense* set of rotations, so any single-qubit unitary can be approximated to arbitrary accuracy.

**What the code does**
* Builds `H·T·H` and verifies `HTH = e^{iπ/8}(cos(π/8)I − i sin(π/8)X)`. Printed matrices **match exactly** . Circuits shown as `─H─T─H─`.
* Builds `T·H·T·H` — note that the code is written *in reverse notation order* (`h; t; h; t` → circuit left-to-right is `H–T–H–T`, which as an operator product is `T·H·T·H`), computes `η = 2·arccos(cos²(π/8)) = 1.09606`, notes `η/π = 0.348886` is irrational (so the rotation angle is never a rational multiple of π → the orbit never closes), and verifies `THTH = Rp(η)` **up to global phase** with the notes' axis `p = (cos π/8, sin π/8, cos π/8)`. Printed **True**.
* **Figure — polar plot "theta vs phi after 1..12 applications of T.H.T.H" :** a 2-D polar scatter/line plot where the radius is the Bloch **θ** and the angle is the Bloch **φ** for 1, 2, …, 12 repetitions of `THTH`, starting from `|+⟩`. It draws a **spiral that keeps stepping around without repeating** — the visual argument for "dense coverage of the sphere ⇒ universality". It is the only polar plot in the notebook.

> **Caveat I verified:** "η/π is irrational" is asserted numerically, not proved; and 12 points can't *show* dense coverage. It's a nice heuristic picture, not a proof (the notes prove it via the algebraic structure of the rotation).

---

## 6. Section 4 — Multi-qubit states and controlled gates 

**Definition under test :** `CU|a,b⟩ = U^{a}|a,b⟩` — apply `U` to the target *only* when the control is 1.

### 6.1 Controlled-H 
The H matrix is wrapped as `UnitaryGate(...).control(1)` and appended as a 2-qubit gate.
* **Circuit figure :** `q0: ─●─` / `q1: ─H─` with the standard filled control dot.
* Printed table for all four inputs (order of the printed vector is Qiskit's little-endian `|00⟩,|01⟩,|10⟩,|11⟩`): `CU|00⟩=|00⟩`, `CU|01⟩=(|01⟩+|11⟩)/√2`, `CU|10⟩=|10⟩`, `CU|11⟩=(|01⟩−|11⟩)/√2`. I re-ran this and it matches exactly — the control is qubit 0, so nothing happens for the `q0=0` inputs.

### 6.2 CNOT / Controlled-X 
* Circuit `q0: ─●─`, `q1: ─⊕─` (drawn twice).
* **Truth table by measurement :** for each `(a,b)` the code prepares the input with `X`s, applies CNOT, measures with 256 shots, and prints the counts:
  `00→{00}`, `01→{10}`, `10→{11}`, `11→{01}` — exactly `b_out = a XOR b`, each a **single 256-count bar** (four histogram figures on "Measurement outcomes"). Doing the truth table *by measured statistics* rather than by matrix algebra is the most important pedagogical upgrade in the notebook.

### 6.3 Bell state 
* Circuit : `q0: ─H─●─`, `q1: ───⊕─`.
* **Figure — state-city plot :** the 3-D "city" of amplitudes, real and imaginary, over the four basis states. Only the `|00⟩` and `|11⟩` skyscrapers of height 1/√2 ≈ 0.707 stand, and they are **blue = purely real, phase 0**. It shows directly that the state is `(|00⟩+|11⟩)/√2` with no relative phase.
* **Figure — Q-sphere :** a globe where each basis state is a node; node size = |amplitude|², node *colour* = phase (colour wheel: 0 → red/pink, π/2 → blue, π → cyan/green, 3π/2 → yellow). You see exactly two equally big nodes, `|00⟩` at the top and `|11⟩` at the bottom, **both the same colour** → the two amplitudes have the same phase. It is the standard "this is a maximally entangled 2-qubit state" picture.
* **Histograms :** "Bell state measurement: only 00 and 11 appear (perfect correlation)", counts `{'00': 2048, '11': 2048}` — the two outcomes that never appear (`01`, `10`) *don't exist at all* in the plot. That absence is the fingerprint of entanglement.

### 6.4 SWAP gate 
* Circuit : `X` on q0 (prepare `|q1q0⟩ = |01⟩`), barrier, then **three CNOTs** `cx(0,1), cx(1,0), cx(0,1)`.
* **Histogram :** counts `{'10': 2048}` — a single spike. Input was q0=1, q1=0; after the swap q0=0, q1=1, which Qiskit reports as the string `'10'`. Perfect exchange.
> **Subtlety worth flagging:** the notes give *three* equivalent three-CNOT identities, one for each pair `(a,b)` of the input. The notebook only implements and verifies the `(|01⟩,|10⟩)` version. Also, the "different" constructions differ only because of the CNOT **direction convention** — don't be confused if the notes' picture looks like a different CNOT ordering.

### 6.5 Control-by-|0⟩ (open control)
* `cx(ctrl_state=0)` draws the **hollow open circle** on the control line and triggers the target when the control is `|0⟩`.
* The code then builds `X–CNOT–X` on the control line (p. 59) and verifies the two operators are **identical** (`True`) — the notes identity.

---

## 7. Section 5 — n-control qubit gates 

### 7.1 Toffoli (CCX) :
* Circuit: `q0: ─●─`, `q1: ─●─`, `q2: ─⊕─`.
* **Truth table by measurement (p. 60):** all 8 inputs, 256 shots each, single spikes: `c_out = c XOR (a AND b)` — the target flips only for `a=b=1`. Seven pages of "Measurement outcomes" histograms each with one bar.

### 7.2 Multi-controlled-X (Cn) :
* `MCXGate(4)` with all four controls preset to `|1⟩` via `X` gates (circuit, the classic 5-line fan-in picture).
* **Histogram :** `{'11111': 2048}` — deterministic flip of the target.
* **Decomposition figure :** `transpile(qc, basis_gates=['u','cx'])` printed as a tall multi-page circuit — the explicit gate-level decomposition of a C⁴X into CNOTs and single-qubit `u` gates (what the notes describe as repeated use).
---

## 8. Section 6 — Controlled-phase and phase kickback 

**Physics :** `CP(α)` adds the phase `e^{iα}` **only** to the `|11⟩` component. Because the phase sits on a *product* state, it can be reinterpreted as a phase on the control → *phase kickback*.

* Circuit : `q0: ─H─●─H─` (control in superposition) and `q1: ─X─P(α)─` (target in `|1⟩`), i.e. `H`, `CP(α)`, `H`.
* **Q-sphere figure (p. 74–76):** two nodes of equal size, `|10⟩` at the equator (control=1… in Qiskit's ordering `|q1q0⟩`) and `|11⟩` at the bottom. Crucially they now have **different colours**: `|10⟩` is blue (phase 0) while `|11⟩` is purple (phase π/3 = 60°) — the phase wheel in the corner lets you read the kicked-back phase α directly off the picture. Printed check: `measured relative phase = 1.0471975511965976 = π/3 = α` ✔.
* **Figure — the phase-kickback interference scan (p. 76–77, plot on p. 101):** α is swept over 25 values from 0 to 2π; for each, the notebook runs `H–CP(α)–H` on the control, measures, and plots `P(control = 1)` vs α. The curve is **1 − cos²(...)**-shaped: it starts at 0, rises to a **maximum of 1.0 at α = π**, and returns to 0 at 2π. This single plot is the whole section: the phase α — which is *completely unobservable* in any Z-basis measurement of the target — becomes a full-swing probability on the control. (This is the mechanism behind phase-estimation and nearly every phase-based quantum algorithm.)
* **Example histogram :** with α = π/2 the outcome is a genuine 50/50 (968 vs 1080). At α = 0 the plot shows a clean 2048–0, and near α = π a clean 0–2048 (the histograms are the 24 individual runs behind the scan).

---

## 9. Section 7 — The no-cloning theorem 

**Physics :** no unitary can copy an arbitrary unknown `|ψ⟩`, because unitaries are linear and the "cloning" map is not.

**The experiment:** the most natural candidate cloner is `CNOT` with the unknown state on the control and a blank `|0⟩` ancilla as target: `U(θ,φ)` then `CX(0,1)`.

* Circuit : `q0: ─U(θ,φ)─●─`, `q1: ───────⊕─`.
* **Fidelity table :** for θ = 0 … π the state fidelity between the circuit output and the ideal `|ψ⟩|ψ⟩` is printed:

| θ | 0 | 0.524 | 1.047 | 1.571 | 2.094 | 2.618 | π |
|---|---|---|---|---|---|---|---|
| fidelity | **1.000** | 0.844 | 0.600 | **0.500** | 0.600 | 0.844 | **1.000** |

  The pattern is exact: **F = (cos³(θ/2) + sin³(θ/2))²**, which equals 1 for the basis states (`|0⟩`, `|1⟩` — classical bits *can* be copied) and drops to its minimum ½ for `|+⟩` (θ = π/2, the maximally "quantum" direction). This is the quantitative form of the theorem, not just a slogan.
* **Histograms (p. 106–109):**
  * "Attempted clone of |0⟩ (basis state): looks like cloning" — a **single 2048 bar**.
  * "Attempted clone of |+⟩ (superposition): entanglement, NOT two copies of |+⟩" — **two ~50/50 bars**. The comment in the notebook is the key: what CNOT actually produces is the *entangled* state `(|00⟩+|11⟩)/√2`; each qubit alone is a maximally mixed 50/50, so no qubit carries "the original `|+⟩`".
* **The decisive test :** apply `H` to *both* qubits (i.e. measure both in the X basis) and look for the "cloning worked" signature, which would be a single bar at `00`. Circuit and histogram show **`{'11': 1000, '00': 1048}`** — a spread. The notebook's own conclusion : *"The histogram shows a spread across outcomes, not a clean spike at 00 — direct experimental confirmation of the no-cloning theorem."*

This is arguably the best-designed section of the notebook: prediction → quantitative curve → histogram → falsification test, with the failure used as the evidence.

---

## 10. Section 8 : Deutsch's algorithm 

**Physics :** with the oracle `U_f|x⟩|y⟩ = |x⟩|y ⊕ f(x)⟩`, the circuit `X(anc) · H ⊗ H · U_f · H(register)` gives `|0⟩` if `f` is constant and `|1⟩` if `f` is balanced — decided with **one** query.

* `deutsch_oracle(kind)` builds all four single-bit oracles: `constant_0` (do nothing), `constant_1` (`X` on the ancilla), `balanced_identity` (`CNOT`), `balanced_not` (`CNOT` then `X`).
* `deutsch_circuit(kind)` wraps the oracle in a labelled `U_f` box with `H`s, barriers and a single measurement on the register qubit.
* **Circuit figures (p. 113–114):** four diagrams, `Deutsch algorithm circuit — f = constant_0 / constant_1 / balanced_identity / balanced_not`, each showing the `Uf[...]` block between two Hadamard layers and the single classical bit `c`. Register qubit `q0` carries the top `H`, the ancilla `q1` carries `X` then `H` (the `|−⟩` preparation of Eq. 4.11), and only `q0` is measured.
  * For `f = constant_0` the oracle box is (correctly) **empty** — "do nothing" *is* the constant-0 oracle, since `f(x) = 0` means `y ⊕ 0 = y`. This is easy to misread as a missing gate; it isn't.
* **Histograms :** four "Measurement outcomes" plots, each a **single 2048 bar**:

| f | counts | verdict |
|---|---|---|
| constant_0 | `{'0': 2048}` | CONSTANT |
| constant_1 | `{'0': 2048}` | CONSTANT |
| balanced_identity | `{'1': 2048}` | BALANCED |
| balanced_not | `{'1': 2048}` | BALANCED |

  Deterministic, correct classification of all four functions **with a single oracle query** — the notebook's clear statement of the quantum advantage (vs. up to 2 classical queries for 1 bit, and N/2+1 in general).

---

## 11. Section 9 : The Quantum Fourier Transform circuit 

### 11.1 Building and drawing the QFT
`qft_manual(n)` implements the notes' construction directly:
```
for j in range(n):
    h(j);                       # Hadamard on j
    for k in range(j+1, n):     # controlled-phase R_k
        cp(2π / 2^(k-j+1), k, j)
swaps to reverse the qubit order
```
* **Circuit figure :** *"Manually-built QFT circuit for n = 3"* (drawn twice) — `H` on the first line, its `P(π/2)`, `P(π/4)` controls from below, then `H`s on the second line, the final `SWAP` between the outer lines. Exactly the textbook QFT ladder. A Qiskit `DeprecationWarning` for `QFT(n)` is also printed here.
* The notebook compares it with `qiskit.circuit.library.QFT(3)`. Printed result (p. 25): **`Manual QFT matches Qiskit library QFT: False`** — this is a *convention* difference (ordering of `H` vs. controlled-phases, and the final SWAP), not a physics error. The library version, drawn on (*"Qiskit built-in QFT (decomposed) for comparison"*, shown twice), has the same inverse-V pattern but accumulates the control-phases in the other direction.

### 11.2 What the QFT does to a basis state
* **Circuit :** two `X` gates to set the register to `|3⟩`, then the manual QFT — *"Prepare basis state, then apply manual QFT"*.
  
* **Figure — amplitude bar chart :** `|amplitude|²` vs basis index y for all 8 outcomes. **All eight bars are exactly equal (0.125)** and the title says it: *"QFT output: equal-magnitude superposition over all y, phase encodes x"*. The picture demonstrates the QFT's defining behaviour: amplitudes are uniform, **all the input information moves into the phases**.
  
* **Q-sphere figure :** eight nodes around the globe of equal size (all |amp|² = 1/8) but **different colours** — the phase wheel shows eight distinct hues, i.e. the phase `2π·x·y/8` winding around.
* **Histogram (p. 30):** *"Measuring QFT output: flat distribution over all 2^n outcomes"*, with circuit *"Full circuit: state prep + QFT + measurement"* shown above it — eight bars of **513, 510, 508, 496, 539, 560, 489, 481** out of 4096 (labels `000`…`111`). Because measurement discards phase, a QFT output always looks flat; the flat histogram is the honest illustration of "phases are invisible in a Z-basis measurement".

### 11.3 Period-finding demo 
* **Text + repeated figures :** *"Period-finding style demonstration (cf. Sec. 4.3 example, Eq. 4.46-4.49). We replicate the worked example of Sec. 4.3.1: a periodic function with period P=2 on n=3 bits..."* — the flat-distribution histogram and full circuit from p. 30 are shown again here before the new code, plus the setup `n_bits = 3; N = 8; P = 2`.
* **Circuit (p. 32):** `H`s on the 3-bit register → a "naive oracle" (`cx(2, ancilla)`) → the manual `QFT` block → measure the register. Title: *"Period-finding style circuit, N = 8, true period P = 2"*.
* **Histogram and conclusion :** counts `{'000': 2084, '001': 2012}` — peaks at the register values 0 and 4 = **N/P**, exactly as the notes' Eq. 4.48 predicts. The printed line repeats the theory: *"Compare with theory (Eq. 4.48): peaks should appear at y = 0 and y = N/2 = 4"*.
> **Read this one carefully:** the plotting labels are Qiskit bit-strings (`'001'`), not integers, so the bar labelled `001` *is* the peak at y = 4 under the notebook's bit-ordering convention. Also, the "oracle" here is not a genuine period-2 oracle — it's a single CNOT chosen to reproduce the worked example's counts, not a general period-finding black box.

---

## 12. Summary table & exercises 

The notebook closes with a **mapping table** (`Notes section → Gate(s) → Notebook section → What's new`), and five **student exercises** :

1. Repeat the phase-kickback scan with a Toffoli-controlled phase (2 controls).
2. Extend Deutsch to 2-bit **Deutsch–Jozsa**; draw and measure all outputs.
3. Change the period-finding demo to P = 4 and confirm peaks move to multiples of N/4.
4. Sweep θ continuously in the no-cloning experiment and plot P(00) in the X basis against the fidelity curve.
5. Build the 4-qubit QFT, draw it, and measure it on a GHZ-like periodic input.

---

## 13. Figure-type index — what you are looking at

| Figure type | Appearance | How to read it |
|---|---|---|
| **Circuit diagram** (`draw('mpl', style='bw')`) | Black-and-white boxes on horizontal wires, `q0/q1/…` labels | Boxes = gates, `●` filled dot = control active on `|1⟩`, `○` open dot = control on `|0⟩`, `⊕` = XOR target, `┼` = SWAP, `▧`=measurement box, `meas` = classical register |
| **Bloch sphere** | Grey globe, magenta arrow, axes x/y/z, `|0⟩`/`|1⟩` poles | Arrow direction = qubit state; rotation about z changes azimuth only; about x tilts it out of the equator; a pure phase gate leaves the arrow pointing at the same spot |
| **Histogram** (`plot_histogram`) | Blue bars, count axis, title | One tall bar = deterministic outcome; two equal bars = 50/50 superposition; **missing** bars = outcomes forbidden by entanglement; flat bars = phase information destroyed by measurement |
| **State-city** (`plot_state_city`) | 3-D skyscraper bars over `\|00⟩…\|11⟩`, two colours = real/imag | Bar height = amplitude, colour = real (blue) vs imaginary; equal heights at `00` and `11` = Bell state |
| **Q-sphere** (`plot_state_qsphere`) | Grey globe with coloured nodes + phase wheel | Node size = probability, node **colour = phase** (wheel: 0 → pink/red, π/2 → blue, π → green, 3π/2 → yellow); equal-phase nodes are the same colour |
| **Polar plot** | Line of points on a circular grid | Radius = Bloch θ, angle = Bloch φ; shows the orbit of repeated `THTH` never closing |
| **Line/scatter scan** | `P(control = 1)` vs α | Rising to 1 at α = π = phase kickback; the full interference curve |
| **Bar chart of probabilities** | `\|amplitude\|²` vs basis index | QFT output: all bars equal → information is in the phases |
| **Tall decomposed circuit** | Very long multi-page circuit | `MCX` expanded into `u` + `cx` by the transpiler |
