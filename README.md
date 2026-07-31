# TorsionLab — Newtonian Three-Body Simulation

TorsionLab generates a daily 1080×1080 H.264 MP4 animation of three mutually gravitating bodies and can publish it to a Telegram channel. Each run creates random planar initial conditions, integrates the equations of motion, renders trajectories, and optionally sends the animation through the Telegram Bot API.

The project is a transparent numerical experiment in **classical Newtonian gravity**, not an orbital prediction service.

## Physical model

The system consists of three **point masses** moving in a two-dimensional plane. Their only interaction is pairwise Newtonian gravitation.

For body $i$, with position $\mathbf r_i$, mass $m_i$, and acceleration $\ddot{\mathbf r}_i$, the code implements:

$$
\ddot{\mathbf r}_i =
G \sum_{j \ne i} m_j
\frac{\mathbf r_j - \mathbf r_i}
{\left|\mathbf r_j - \mathbf r_i\right|^3}.
$$

Equivalently, the force exerted on body $i$ by body $j$ is:

$$
\mathbf F_{ij} =
G\frac{m_i m_j}{\left|\mathbf r_j-\mathbf r_i\right|^3}
(\mathbf r_j-\mathbf r_i).
$$

Newton's second law is:

$$
m_i \frac{d^2 \mathbf r_i}{dt^2} = \sum_j \mathbf F_{ij}
$$

This relation produces the acceleration equation. The project uses normalized units with $G=1$, so the animation does not claim SI units such as kilograms or metres.

### Newton's laws in the model

1. **Newton's second law:** net gravitational force determines acceleration.
2. **Universal gravitation:** each pair attracts with an inverse-square force; the vector acceleration has an inverse-cube denominator.
3. **Newton's third law:** pair forces are equal and opposite. Initial velocities are adjusted so total initial linear momentum is zero.

## Initial conditions

Every generated system contains exactly three bodies.

- Each mass is sampled independently from the uniform interval $[0.6, 1.8]$.
- Each initial position is sampled uniformly from $[-5,5]\times[-5,5]$.
- A configuration is accepted only when every initial pair separation is at least $1.2$.
- Each raw velocity component is drawn from a normal distribution with standard deviation $0.45$.
- The mass-weighted mean velocity is removed:

$$
\mathbf v_i \leftarrow \mathbf v_i -
\frac{\sum_k m_k\mathbf v_k}{\sum_k m_k}.
$$

Therefore the initial total linear momentum

$$
\mathbf P = \sum_i m_i\mathbf v_i
$$

is zero up to floating-point arithmetic. Positions are deliberately **not** shifted to the centre of mass, so the bodies begin at visibly different random locations in the GIF.

The Telegram caption reports masses and initial momenta:

$$
\mathbf p_i = m_i\mathbf v_i.
$$

## Numerical integration

The equations are solved with **velocity-Verlet**, a time-reversible symplectic method suited to conservative mechanical systems.

For timestep $\Delta t$:

$$
\mathbf r_{n+1} =
\mathbf r_n + \mathbf v_n\Delta t +
\frac{1}{2}\mathbf a_n\Delta t^2,
$$

$$
\mathbf a_{n+1} = \mathbf a(\mathbf r_{n+1}),
$$

$$
\mathbf v_{n+1} =
\mathbf v_n + \frac{1}{2}
(\mathbf a_n+\mathbf a_{n+1})\Delta t.
$$

| Parameter | Value | Meaning |
| --- | ---: | --- |
| Bodies | 3 | Point masses in a plane |
| $G$ | 1 | Normalized gravitational constant |
| Duration | 30 | Simulation time units |
| $\Delta t$ | 0.002 | Integrator timestep |
| Captured frames | 360 | Approximate MP4 sampling target |
| Trail length | 100 frames | Recent trajectory drawn for each body |

The integrator runs at the small timestep; only selected states are stored for rendering. The renderer does not alter the physical calculation.

## Conservation diagnostics

For every published three-body run, the program evaluates the saved states and reports two numerical diagnostics in the Telegram caption:

- **Relative energy drift:** the largest sampled value of $|E(t)-E(0)|/|E(0)|$.
- **Angular-momentum drift:** the largest sampled value of $|L_z(t)-L_z(0)|$.

The total mechanical energy is calculated as

$$
E = \frac{1}{2}\sum_i m_i |\mathbf v_i|^2
    - G\sum_{i<j}\frac{m_i m_j}{|\mathbf r_j-\mathbf r_i|},
$$

and the planar angular momentum is

$$
L_z = \sum_i m_i (x_i v_{y,i} - y_i v_{x,i}).
$$

These values are diagnostics of the **numerical integration**, not corrections applied to the simulation. Small drift is expected from finite timesteps; a large drift indicates that the timestep or initial conditions should be investigated.

Each publication also records its random **seed**, integration method, timestep, and duration. Re-run the simulation by passing that seed to `random_initial_conditions(seed=...)`.

## Two-body validation

The script includes a known circular orbit for two equal point masses. It uses a separation of one normalized distance unit and the corresponding circular speed $v=\sqrt{1/2}$. Run the validation without generating or posting a GIF:

~~~
python main.py --validate
~~~

The command prints the maximum relative energy drift and absolute angular-momentum drift over approximately one orbit. This is a sanity check for the velocity-Verlet implementation; it is not a proof of accuracy for every chaotic three-body configuration.

## Close encounters and limitations

The point-mass force is singular at zero separation. This project does **not** soften the force, merge bodies, or model collisions.

If any pair comes closer than

$$
r_\mathrm{min}=0.02,
$$

integration stops with a **CloseEncounter** event. That random system is discarded, and the program tries another one, up to 20 attempts.

The model does not include:

- collisions, finite radii, tidal deformation, or mass transfer;
- general relativity or gravitational radiation;
- external fields, drag, or relativistic corrections;
- three-dimensional motion;
- adaptive timesteps or long-term error analysis.

The MP4 animation is therefore a qualitative numerical visualization under the stated assumptions.

## Program flow

~~~
random initial conditions
        ↓
Newtonian accelerations
        ↓
velocity-Verlet integration
        ↓
MP4 rendering with positions and recent trails
        ↓
Telegram sendAnimation request (when configured)
~~~

The active entry point is **compact_main()** in **main.py**. It generates a valid system, saves **daily_three_body.mp4**, builds a scientific caption with masses, momenta, seed, and conservation diagnostics, and calls the Telegram publisher.

## Installation

Requirements:

- Python 3.11 or newer
- Packages in **requirements.txt**

Install dependencies:

~~~
python -m pip install -r requirements.txt
~~~

Run locally:

~~~
python main.py
~~~

The MP4 video is written beside the script as **daily_three_body.mp4**.

## Telegram publishing

Posting occurs only when both environment variables exist:

~~~
TELEGRAM_BOT_TOKEN
TELEGRAM_CHANNEL_ID
~~~

The bot must be a channel administrator with permission to post messages. The code calls Telegram's **sendAnimation** endpoint and attaches the H.264 MP4 as the animation file. The caption uses Telegram HTML formatting for readable scientific metadata.

Never place a token in source code, commit it to GitHub, or share it in chat. Revoke an exposed token through @BotFather and issue a replacement.

## Daily automation on Windows

Create a daily Windows Task Scheduler task that runs:

~~~
python main.py
~~~

Set the project directory as the working directory. The process runs only while it generates and posts, then exits; it is not a permanently active background service.

## Repository structure

| File | Purpose |
| --- | --- |
| **main.py** | Simulation, integration, rendering, and Telegram publishing |
| **test_physics.py** | Regression test for the two-body validation orbit |
| **requirements.txt** | Python dependencies |
| **README.md** | Scientific and operational documentation |

## License

No license has been selected yet. Add one before accepting external contributions or defining reuse terms.
