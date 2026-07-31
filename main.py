# =============================================================================
# TorsionLab: quick reading guide
#
# This file has four responsibilities:
#   1. Create a random Newtonian three-body system.
#   2. Integrate its motion with the velocity-Verlet method.
#   3. Render the sampled positions as a GIF animation.
#   4. Optionally post that GIF to a Telegram channel.
#
# Suggested order for a first read:
#   Simulation              — groups the input data for one experiment.
#   random_initial_conditions — chooses masses, positions, and velocities.
#   accelerations           — implements Newtonian gravity.
#   solve                   — advances the system through time.
#   compact_render          — creates the Telegram-friendly GIF.
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

import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import requests
from matplotlib.animation import FuncAnimation, PillowWriter


G = 1.0
COLLISION_TOLERANCE = 0.02


class CloseEncounter(Exception):
    """Raised before the point-mass model reaches its singularity."""


@dataclass
class Simulation:
    masses: np.ndarray
    positions: np.ndarray
    velocities: np.ndarray
    duration: float = 30.0
    time_step: float = 0.002


def random_initial_conditions(seed: Optional[int] = None) -> Simulation:
    """Create a bounded, visually useful initial configuration."""
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
    return Simulation(masses, positions, velocities)


def accelerations(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    acceleration = np.zeros_like(positions)
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            delta = positions[j] - positions[i]
            distance_sq = float(delta @ delta)
            if distance_sq < COLLISION_TOLERANCE**2:
                raise CloseEncounter("Тіла досягли порогу близького зближення.")
            acceleration[i] += G * masses[j] * delta / distance_sq ** 1.5
    return acceleration


def solve(simulation: Simulation, frames: int = 360) -> np.ndarray:
    """Integrate with velocity-Verlet (a stable symplectic method)."""
    steps = round(simulation.duration / simulation.time_step)
    capture_every = max(1, steps // frames)
    positions = simulation.positions.copy()
    velocities = simulation.velocities.copy()
    acceleration = accelerations(positions, simulation.masses)
    samples = [positions.copy()]

    for step in range(steps):
        positions += velocities * simulation.time_step + 0.5 * acceleration * simulation.time_step**2
        new_acceleration = accelerations(positions, simulation.masses)
        velocities += 0.5 * (acceleration + new_acceleration) * simulation.time_step
        acceleration = new_acceleration
        if step % capture_every == 0:
            samples.append(positions.copy())
    return np.asarray(samples)


def render(trajectory: np.ndarray, simulation: Simulation, destination: Path) -> Path:
    """Render a research-style animation with the governing equation and data."""
    masses = simulation.masses
    colors = ("#e76f51", "#2a9d8f", "#457b9d")
    limit = max(2.0, float(np.abs(trajectory).max()) * 1.15)
    fig = plt.figure(figsize=(8, 8), facecolor="#f8f9fa")
    grid = fig.add_gridspec(2, 1, height_ratios=(6.3, 1.3), hspace=0.28,
                            top=0.91, bottom=0.055)
    ax = fig.add_subplot(grid[0])
    details = fig.add_subplot(grid[1])
    ax.set_facecolor("#ffffff")
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal")
    ax.set_xlabel("x  [нормовані одиниці]")
    ax.set_ylabel("y  [нормовані одиниці]")
    ax.grid(True, color="#d7dde3", linewidth=0.7)
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
        response = requests.post(url, data={"chat_id": channel, "caption": caption}, files={"animation": animation}, timeout=60)
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
    output = Path(__file__).resolve().parent / "daily_three_body.gif"
    render(trajectory, simulation, output)
    mass_text = ", ".join(f"{mass:.2f}" for mass in simulation.masses)
    caption = (f"Ньютонівська задача трьох тіл — {datetime.now():%d.%m.%Y}\n"
               f"Маси тіл: {mass_text}\nМодель: точкові маси, без зіткнень і згладжування сили.")
    publish_to_telegram(output, caption)
    print(f"Done: {output.resolve()}")


def compact_render(trajectory: np.ndarray, simulation: Simulation, destination: Path) -> Path:
    """Render a clean, Telegram-friendly orbit visualization without text panels."""
    colors = ("#e76f51", "#2a9d8f", "#457b9d")
    limit = max(2.0, float(np.abs(trajectory).max()) * 1.15)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=140, facecolor="#f8f9fa")
    fig.subplots_adjust(left=0.12, right=0.96, top=0.965, bottom=0.11)
    ax.set_facecolor("#ffffff")
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal")
    ax.set_xlabel("x  [нормовані одиниці]")
    ax.set_ylabel("y  [нормовані одиниці]")
    ax.grid(True, color="#d7dde3", linewidth=0.7)
    for spine in ax.spines.values():
        spine.set_color("#495057")
    points = [ax.plot([], [], "o", color=color, ms=7 + 5 * mass)[0]
              for color, mass in zip(colors, simulation.masses)]
    trails = [ax.plot([], [], "-", color=color, lw=1.25, alpha=0.75)[0] for color in colors]
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


def build_caption(simulation: Simulation) -> str:
    masses = "  ".join(f"m{index}={mass:.3f}" for index, mass in enumerate(simulation.masses, start=1))
    momenta = simulation.masses[:, np.newaxis] * simulation.velocities
    momentum_text = "  ".join(
        f"p{index}=({momentum[0]:+.3f}, {momentum[1]:+.3f})"
        for index, momentum in enumerate(momenta, start=1)
    )
    return (f"Torsivane Lab — {datetime.now():%d.%m.%Y}\n"
            f"Маси: {masses}\n"
            f"Початкові імпульси: {momentum_text}")


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
    output = Path(__file__).resolve().parent / "daily_three_body.gif"
    compact_render(trajectory, simulation, output)
    publish_to_telegram(output, build_caption(simulation))
    print(f"Done: {output.resolve()}")


if __name__ == "__main__":
    compact_main()
