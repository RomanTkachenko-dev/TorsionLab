"""Small regression test for the numerical integrator.

Run with: python -m unittest test_physics.py
"""

import unittest

from main import validate_two_body_orbit


class TwoBodyValidationTests(unittest.TestCase):
    """A circular two-body orbit should preserve its invariants very closely."""

    def test_velocity_verlet_conserves_invariants_for_one_orbit(self) -> None:
        report = validate_two_body_orbit()
        self.assertLess(report.energy_relative_drift, 1e-8)
        self.assertLess(report.angular_momentum_absolute_drift, 1e-10)


if __name__ == "__main__":
    unittest.main()
