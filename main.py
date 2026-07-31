# =============================================================================
# TorsionLab: quick reading guide
#
# This file has four responsibilities:
#   1. Create a random Newtonian three-body system.
#   2. Integrate its motion with the velocity-Verlet method.
#   3. Render the sampled positions as a 1080×1080 H.264 MP4 animation.
#   4. Optionally post that animation to a Telegram channel.
#
# Suggested order for a first read:
#   Simulation              — groups the input data for one experiment.
#   random_initial_conditions — chooses masses, positions, and velocities.
#   accelerations           — implements Newtonian gravity.
#   solve                   — advances the system through time.
#   compact_render          — creates the Telegram-friendly MP4 video.
#   build_caption           — formats masses and initial momenta for Telegram.
#   compact_main            — the active high-level workflow.
#
# Important conventions:
#   * G = 1 uses normalized units rather than SI units.
#   * Bodies are point masses; collisions and mergers are not modelled.
#   * A close encounter aborts that random attempt before the force singularity.
#   * Lines ending in `_main` coordinate work; the active entry point is
#     the compact_main() call at the bottom of the file.
# =============================================================================

"""Generate a daily three-body animation and optionally post it to Telegram."""

from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import requests
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter


G = 1.0
COLLISION_TOLERANCE = 0.02


def application_directory() -> Path:
    """Return the folder containing the script during development or the EXE after packaging."""
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def configure_bundled_ffmpeg() -> None:
    """Tell Matplotlib where to find FFmpeg when this file runs as a packaged EXE."""
    bundle_directory = Path(getattr(sys, "_MEIPASS", application_directory()))
    bundled_ffmpeg = bundle_directory / "ffmpeg.exe"
    if bundled_ffmpeg.is_file():
        plt.rcParams["animation.ffmpeg_path"] = str(bundled_ffmpeg)


configure_bundled_ffmpeg()


class CloseEncounter(Exception):
    """Raised before the point-mass model reaches its singularity."""


@dataclass
class Simulation:
    masses: np.ndarray
    positions: np.ndarray
    velocities: np.ndarray
    duration: float = 30.0
    time_step: float = 0.002
    seed: Optional[int] = None


@dataclass
class Solution:
    """States sampled from the integrator for rendering and diagnostics."""

    positions: np.ndarray
    velocities: np.ndarray


@dataclass
class Diagnostics:
    """Conservation metrics evaluated at the saved integration states."""

    energy_initial: float
    energy_relative_drift: float
    angular_momentum_initial: float
    angular_momentum_absolute_drift: float


@dataclass
class OutcomeAnalysis:
    """Physical summary used for the research log and future statistics."""

    closest_pair: str
    minimum_separation: float
    closest_approach_time: float
    final_bound_pair: str
    final_pair_relative_energy: float
    escape_candidate: str
    escape_specific_energy: float
    escape_radial_speed: float


def random_initial_conditions(seed: Optional[int] = None) -> Simulation:
    """Create a bounded, visually useful initial configuration."""
    # Keep the actual seed with the simulation so a published run is reproducible.
    if seed is None:
        seed = int.from_bytes(os.urandom(8), "big")
    rng = np.random.default_rng(seed)
    masses = rng.uniform(0.6, 1.8, size=3)
    while True:
        positions = rng.uniform(-5.0, 5.0, size=(3, 2))
        separations = [np.linalg.norm(positions[i] - positions[j])
                       for i in range(3) for j in range(i + 1, 3)]
        if min(separations) >= 1.2:
            break
    velocities = rng.normal(0, 0.45, size=(3, 2))
    velocities -= np.average(velocities, axis=0, weights=masses)
    return Simulation(masses, positions, velocities, seed=seed)


def accelerations(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Return Newtonian gravitational acceleration for every body."""
    acceleration = np.zeros_like(positions)
    body_count = len(masses)
    for i in range(body_count):
        for j in range(body_count):
            if i == j:
                continue
            delta = positions[j] - positions[i]
            distance_sq = float(delta @ delta)
            if distance_sq < COLLISION_TOLERANCE**2:
                raise CloseEncounter("Тіла досягли порогу близького зближення.")
            acceleration[i] += G * masses[j] * delta / distance_sq ** 1.5
    return acceleration


def solve(simulation: Simulation, frames: int = 360) -> Solution:
    """Integrate with velocity-Verlet (a stable symplectic method)."""
    steps = round(simulation.duration / simulation.time_step)
    capture_every = max(1, steps // frames)
    positions = simulation.positions.copy()
    velocities = simulation.velocities.copy()
    acceleration = accelerations(positions, simulation.masses)
    position_samples = [positions.copy()]
    velocity_samples = [velocities.copy()]

    for step in range(steps):
        positions += velocities * simulation.time_step + 0.5 * acceleration * simulation.time_step**2
        new_acceleration = accelerations(positions, simulation.masses)
        velocities += 0.5 * (acceleration + new_acceleration) * simulation.time_step
        acceleration = new_acceleration
        if step % capture_every == 0:
            position_samples.append(positions.copy())
            velocity_samples.append(velocities.copy())
    return Solution(np.asarray(position_samples), np.asarray(velocity_samples))


def diagnostics(solution: Solution, simulation: Simulation) -> Diagnostics:
    """Measure conservation of total energy and z-angular momentum."""
    positions = solution.positions
    velocities = solution.velocities
    masses = simulation.masses

    kinetic = 0.5 * np.sum(masses[np.newaxis, :, np.newaxis] * velocities**2, axis=(1, 2))
    potential = np.zeros(len(positions))
    for i in range(len(masses)):
        for j in range(i + 1, len(masses)):
            separation = np.linalg.norm(positions[:, j] - positions[:, i], axis=1)
            potential -= G * masses[i] * masses[j] / separation
    total_energy = kinetic + potential

    angular_momentum = np.sum(
        masses[np.newaxis, :] * (positions[:, :, 0] * velocities[:, :, 1]
                                 - positions[:, :, 1] * velocities[:, :, 0]),
        axis=1,
    )

    energy_scale = max(abs(total_energy[0]), np.finfo(float).tiny)
    return Diagnostics(
        energy_initial=float(total_energy[0]),
        energy_relative_drift=float(np.max(np.abs(total_energy - total_energy[0])) / energy_scale),
        angular_momentum_initial=float(angular_momentum[0]),
        angular_momentum_absolute_drift=float(np.max(np.abs(angular_momentum - angular_momentum[0]))),
    )


def analyse_outcome(solution: Solution, simulation: Simulation) -> OutcomeAnalysis:
    """Extract simple, reproducible descriptors from one completed run.

    A negative relative two-body energy at the final saved state marks a
    *bound-pair candidate*. An *escape candidate* is the remaining body when
    it is moving away from that pair, has positive pair-relative specific
    energy, and lies farther from the pair's centre of mass than the pair is
    from itself. These are finite-window descriptors, not permanent outcomes.
    """
    positions = solution.positions
    velocities = solution.velocities
    masses = simulation.masses
    pairs = [(i, j) for i in range(len(masses)) for j in range(i + 1, len(masses))]

    closest_pair = ""
    minimum_separation = np.inf
    closest_frame = 0
    relative_energies: list[tuple[float, tuple[int, int]]] = []
    for i, j in pairs:
        separations = np.linalg.norm(positions[:, j] - positions[:, i], axis=1)
        frame = int(np.argmin(separations))
        if separations[frame] < minimum_separation:
            minimum_separation = float(separations[frame])
            closest_pair = f"m{i + 1}-m{j + 1}"
            closest_frame = frame

        reduced_mass = masses[i] * masses[j] / (masses[i] + masses[j])
        relative_speed_sq = float(np.sum((velocities[-1, i] - velocities[-1, j]) ** 2))
        final_distance = float(np.linalg.norm(positions[-1, j] - positions[-1, i]))
        relative_energy = 0.5 * reduced_mass * relative_speed_sq - G * masses[i] * masses[j] / final_distance
        relative_energies.append((relative_energy, (i, j)))

    final_energy, (i, j) = min(relative_energies, key=lambda item: item[0])
    final_pair = f"m{i + 1}-m{j + 1}" if final_energy < 0 else "none"
    escape_candidate = "none"
    escape_specific_energy = float("nan")
    escape_radial_speed = float("nan")

    if final_energy < 0:
        escaping_index = next(index for index in range(len(masses))
                              if index not in (i, j))
        pair_mass = masses[i] + masses[j]
        pair_position = (masses[i] * positions[-1, i] + masses[j] * positions[-1, j]) / pair_mass
        pair_velocity = (masses[i] * velocities[-1, i] + masses[j] * velocities[-1, j]) / pair_mass
        separation_vector = positions[-1, escaping_index] - pair_position
        separation = float(np.linalg.norm(separation_vector))
        relative_velocity = velocities[-1, escaping_index] - pair_velocity
        escape_radial_speed = float(separation_vector @ relative_velocity / separation)
        escape_specific_energy = (
            0.5 * float(relative_velocity @ relative_velocity)
            - G * (pair_mass + masses[escaping_index]) / separation
        )
        pair_separation = float(np.linalg.norm(positions[-1, j] - positions[-1, i]))
        if (escape_specific_energy > 0 and escape_radial_speed > 0
                and separation > pair_separation):
            escape_candidate = f"m{escaping_index + 1}"
    time_per_frame = simulation.duration / (len(positions) - 1)
    return OutcomeAnalysis(
        closest_pair=closest_pair,
        minimum_separation=minimum_separation,
        closest_approach_time=closest_frame * time_per_frame,
        final_bound_pair=final_pair,
        final_pair_relative_energy=float(final_energy),
        escape_candidate=escape_candidate,
        escape_specific_energy=escape_specific_energy,
        escape_radial_speed=escape_radial_speed,
    )


def record_simulation(simulation: Simulation, report: Diagnostics, outcome: OutcomeAnalysis) -> Path:
    """Append one reproducible observation to the local long-term CSV dataset."""
    destination = application_directory() / "simulation_records.csv"
    momenta = simulation.masses[:, np.newaxis] * simulation.velocities
    row = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seed": simulation.seed,
        "method": "velocity-Verlet",
        "duration": simulation.duration,
        "time_step": simulation.time_step,
        "m1": simulation.masses[0], "m2": simulation.masses[1], "m3": simulation.masses[2],
        "r1_x": simulation.positions[0, 0], "r1_y": simulation.positions[0, 1],
        "r2_x": simulation.positions[1, 0], "r2_y": simulation.positions[1, 1],
        "r3_x": simulation.positions[2, 0], "r3_y": simulation.positions[2, 1],
        "p1_x": momenta[0, 0], "p1_y": momenta[0, 1],
        "p2_x": momenta[1, 0], "p2_y": momenta[1, 1],
        "p3_x": momenta[2, 0], "p3_y": momenta[2, 1],
        "closest_pair": outcome.closest_pair,
        "minimum_separation": outcome.minimum_separation,
        "closest_approach_time": outcome.closest_approach_time,
        "final_bound_pair_candidate": outcome.final_bound_pair,
        "final_pair_relative_energy": outcome.final_pair_relative_energy,
        "escape_candidate": outcome.escape_candidate,
        "escape_specific_energy": outcome.escape_specific_energy,
        "escape_radial_speed": outcome.escape_radial_speed,
        "energy_relative_drift": report.energy_relative_drift,
        "angular_momentum_absolute_drift": report.angular_momentum_absolute_drift,
    }
    fieldnames = list(row)
    if destination.exists():
        with destination.open("r", newline="", encoding="utf-8") as file:
            reader = csv.DictReader(file)
            existing_rows = list(reader)
            existing_fieldnames = reader.fieldnames or []
        if existing_fieldnames != fieldnames:
            with destination.open("w", newline="", encoding="utf-8") as file:
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                for existing_row in existing_rows:
                    writer.writerow({name: existing_row.get(name, "") for name in fieldnames})
                writer.writerow(row)
            return destination
    with destination.open("a", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        if not destination.exists() or destination.stat().st_size == 0:
            writer.writeheader()
        writer.writerow(row)
    return destination


def render(solution: Solution, simulation: Simulation, destination: Path) -> Path:
    """Render a research-style animation with the governing equation and data."""
    masses = simulation.masses
    colors = ("#e76f51", "#2a9d8f", "#457b9d")
    trajectory = solution.positions
    limit = max(2.0, float(np.abs(trajectory).max()) * 1.15)
    fig = plt.figure(figsize=(8, 8), facecolor="#f8f9fa")
    grid = fig.add_gridspec(2, 1, height_ratios=(6.3, 1.3), hspace=0.28,
                            top=0.91, bottom=0.055)
    ax = fig.add_subplot(grid[0])
    details = fig.add_subplot(grid[1])
    ax.set_facecolor("#ffffff")
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(True, color="#aeb7c2", linewidth=0.9, alpha=0.85)
    fig.suptitle("Ньютонівська задача трьох тіл", y=0.975, fontsize=15, fontweight="bold")
    fig.text(0.5, 0.942, "чисельне інтегрування: velocity-Verlet", ha="center",
             fontsize=8.5, color="#495057")
    points = [ax.plot([], [], "o", color=c, ms=7 + 5 * m)[0] for c, m in zip(colors, masses)]
    trails = [ax.plot([], [], "-", color=c, lw=1, alpha=0.7)[0] for c in colors]

    details.axis("off")
    equation = (r"$\ddot{\mathbf{r}}_i = G \sum_{j \ne i} m_j "
                r"\frac{\mathbf{r}_j-\mathbf{r}_i}{|\mathbf{r}_j-\mathbf{r}_i|^3}$")
    details.text(0.5, 0.77, equation, ha="center", va="center", fontsize=10)
    mass_text = "   ".join(f"m{index}={mass:.3f}" for index, mass in enumerate(masses, start=1))
    initial_momenta = simulation.masses[:, np.newaxis] * simulation.velocities
    momentum_text = "   ".join(
        f"p{index}=({momentum[0]:+.3f}, {momentum[1]:+.3f})"
        for index, momentum in enumerate(initial_momenta, start=1)
    )
    details.text(0.5, 0.43, f"G=1     {mass_text}     Δt={simulation.time_step:.4f}",
                 ha="center", va="center", fontsize=8.2, family="monospace", color="#343a40")
    details.text(0.5, 0.12, f"Початковий імпульс:     {momentum_text}",
                 ha="center", va="center", fontsize=7.7, family="monospace", color="#343a40")
    time_label = ax.text(0.02, 0.03, "", transform=ax.transAxes, fontsize=10,
                         bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#adb5bd"})

    def update(frame: int):
        start = max(0, frame - 100)
        for body in range(3):
            points[body].set_data([trajectory[frame, body, 0]], [trajectory[frame, body, 1]])
            trails[body].set_data(trajectory[start:frame + 1, body, 0], trajectory[start:frame + 1, body, 1])
        time_label.set_text(f"t = {frame * simulation.duration / (len(trajectory) - 1):.2f}")
        return [*points, *trails, time_label]

    animation = FuncAnimation(fig, update, frames=len(trajectory), interval=33, blit=True)
    animation.save(destination, writer=PillowWriter(fps=30))
    plt.close(fig)
    return destination


def publish_to_telegram(file_path: Path, caption: str) -> None:
    token = os.environ.get("TELEGRAM_BOT_TOKEN")
    channel = os.environ.get("TELEGRAM_CHANNEL_ID")
    if not token or not channel:
        print("Telegram credentials are absent; animation was generated locally.")
        return
    url = f"https://api.telegram.org/bot{token}/sendAnimation"
    with file_path.open("rb") as animation:
        response = requests.post(
            url,
            data={"chat_id": channel, "caption": caption, "parse_mode": "HTML"},
            files={"animation": animation},
            timeout=60,
        )
    response.raise_for_status()


def main() -> None:
    for _ in range(20):
        simulation = random_initial_conditions()
        try:
            trajectory = solve(simulation)
            break
        except CloseEncounter:
            continue
    else:
        raise RuntimeError("Не вдалося згенерувати систему без близького зближення.")
    output = application_directory() / "daily_three_body.gif"
    render(trajectory, simulation, output)
    mass_text = ", ".join(f"{mass:.2f}" for mass in simulation.masses)
    caption = (f"Ньютонівська задача трьох тіл — {datetime.now():%d.%m.%Y}\n"
               f"Маси тіл: {mass_text}\nМодель: точкові маси, без зіткнень і згладжування сили.")
    publish_to_telegram(output, caption)
    print(f"Done: {output.resolve()}")


def compact_render(solution: Solution, simulation: Simulation, destination: Path) -> Path:
    """Render a 1080×1080 H.264 video for sharp Telegram playback."""
    trajectory = solution.positions
    colors = ("#e76f51", "#2a9d8f", "#457b9d")
    limit = max(2.0, float(np.abs(trajectory).max()) * 1.15)
    # 8 in × 135 dpi = 1080 px on each side.
    fig, ax = plt.subplots(figsize=(8, 8), dpi=135, facecolor="#f8f9fa")
    fig.subplots_adjust(left=0.12, right=0.96, top=0.965, bottom=0.11)
    ax.set_facecolor("#ffffff")
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(True, color="#aeb7c2", linewidth=0.9, alpha=0.85)
    for spine in ax.spines.values():
        spine.set_color("#495057")
    # Crosses remain fixed at the launch points.  The two trail layers show the
    # complete history faintly while keeping the most recent motion prominent.
    start_markers = [ax.plot(
        [trajectory[0, body, 0]], [trajectory[0, body, 1]],
        marker="x", color=color, ms=8, mew=1.5, alpha=0.95, linestyle="None",
    )[0] for body, color in enumerate(colors)]
    full_trails = [ax.plot([], [], "-", color=color, lw=0.9, alpha=0.28)[0] for color in colors]
    recent_trails = [ax.plot([], [], "-", color=color, lw=1.7, alpha=0.95)[0] for color in colors]
    points = [ax.plot([], [], "o", color=color, ms=7 + 5 * mass)[0]
              for color, mass in zip(colors, simulation.masses)]
    time_label = ax.text(0.02, 0.03, "", transform=ax.transAxes, fontsize=10,
                         bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#adb5bd"})

    def update(frame: int):
        start = max(0, frame - 100)
        for body in range(3):
            points[body].set_data([trajectory[frame, body, 0]], [trajectory[frame, body, 1]])
            full_trails[body].set_data(trajectory[:frame + 1, body, 0], trajectory[:frame + 1, body, 1])
            recent_trails[body].set_data(trajectory[start:frame + 1, body, 0], trajectory[start:frame + 1, body, 1])
        time_label.set_text(f"t = {frame * simulation.duration / (len(trajectory) - 1):.2f}")
        return [*start_markers, *full_trails, *recent_trails, *points, time_label]

    animation = FuncAnimation(fig, update, frames=len(trajectory), interval=33, blit=True)
    writer = FFMpegWriter(
        fps=30,
        codec="libx264",
        bitrate=-1,
        extra_args=["-crf", "18", "-preset", "medium", "-pix_fmt", "yuv420p"],
    )
    animation.save(destination, writer=writer)
    plt.close(fig)
    return destination


def build_caption(simulation: Simulation, report: Diagnostics, outcome: OutcomeAnalysis) -> str:
    """Build a reproducible, compact scientific record for the Telegram post."""
    momenta = simulation.masses[:, np.newaxis] * simulation.velocities
    mass_text = "  ".join(
        f"m{index} = {mass:.3f}" for index, mass in enumerate(simulation.masses, start=1)
    )
    momentum_lines = "\n".join(
        f"p{index} = ({momentum[0]:+.3f}, {momentum[1]:+.3f})"
        for index, momentum in enumerate(momenta, start=1)
    )
    return (f"<b>Torsivane Lab</b> · {datetime.now():%d.%m.%Y}\n"
            f"<pre>{mass_text}\n"
            f"{momentum_lines}\n\n"
            f"integration\n"
            f"  seed = {simulation.seed}\n"
            f"  method       = velocity-Verlet\n"
            f"  step Δt      = {simulation.time_step:g}\n"
            f"  window T     = {simulation.duration:g}\n"
            f"outcome (within T)\n"
            f"  closest pass = {outcome.closest_pair}, r = {outcome.minimum_separation:.3f}, t = {outcome.closest_approach_time:.2f}\n"
            f"  bound pair   = {outcome.final_bound_pair}\n"
            f"  escape cand. = {outcome.escape_candidate}\n"
            f"numerical check\n"
            f"  max |ΔE/E₀| = {report.energy_relative_drift:.2e}\n"
            f"  max |ΔL|    = {report.angular_momentum_absolute_drift:.2e}</pre>")


def validate_two_body_orbit() -> Diagnostics:
    """Check a known circular two-body orbit and return conservation diagnostics."""
    orbital_speed = np.sqrt(0.5)
    simulation = Simulation(
        masses=np.array([1.0, 1.0]),
        positions=np.array([[-0.5, 0.0], [0.5, 0.0]]),
        velocities=np.array([[0.0, -orbital_speed], [0.0, orbital_speed]]),
        duration=2 * np.pi / np.sqrt(2),
        time_step=0.002,
    )
    return diagnostics(solve(simulation), simulation)


def compact_main() -> None:
    for _ in range(20):
        simulation = random_initial_conditions()
        try:
            trajectory = solve(simulation)
            break
        except CloseEncounter:
            continue
    else:
        raise RuntimeError("Unable to generate a system without a close encounter.")
    output = application_directory() / "daily_three_body.mp4"
    compact_render(trajectory, simulation, output)
    report = diagnostics(trajectory, simulation)
    outcome = analyse_outcome(trajectory, simulation)
    record_path = record_simulation(simulation, report, outcome)
    publish_to_telegram(output, build_caption(simulation, report, outcome))
    print(f"Done: {output.resolve()}")
    print(f"Recorded: {record_path.resolve()}")
    print(f"max relative energy drift = {report.energy_relative_drift:.2e}")
    print(f"max angular-momentum drift = {report.angular_momentum_absolute_drift:.2e}")


if __name__ == "__main__":
    if "--validate" in sys.argv:
        report = validate_two_body_orbit()
        print("Two-body circular-orbit validation")
        print(f"max relative energy drift = {report.energy_relative_drift:.2e}")
        print(f"max angular-momentum drift = {report.angular_momentum_absolute_drift:.2e}")
    else:
        compact_main()
