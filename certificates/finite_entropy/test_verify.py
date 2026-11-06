import copy
import json
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
import verify as v


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case = v.load_case()

    def test_published_enclosure(self):
        r = v.certify(self.case, enabled=True)
        self.assertGreaterEqual(F(r['lower']), F('4.9836437909'))
        self.assertLessEqual(F(r['upper']), F('4.9836437910'))
        self.assertGreaterEqual(F(r['dual_input_gap_upper']), 0)
        self.assertGreaterEqual(F(r['dual_reset_gap_upper']), 0)

    def test_requires_opt_in(self):
        with self.assertRaises(ValueError):
            v.certify(self.case)

    def test_bad_mass(self):
        c = copy.deepcopy(self.case)
        c['p'][0] = str(int(c['p'][0]) + 1)
        with self.assertRaises(ValueError):
            v.validate(c)

    def test_bad_mean_with_mass_preserved(self):
        c = copy.deepcopy(self.case)
        c['p'][0] = str(int(c['p'][0]) - 1)
        c['p'][1] = str(int(c['p'][1]) + 1)
        with self.assertRaises(ValueError):
            v.validate(c)

    def test_zero_negative_and_float_rejected(self):
        for bad in ['0', '-1', 1.25, True, '1.2', 'nan']:
            c = copy.deepcopy(self.case)
            c['t'][2] = bad
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                v.validate(c)

    def test_parameters_not_silently_changed(self):
        for key, value in [('E', 8.0), ('dimension', 320), ('lambda', '0.6'), ('c', '1/10')]:
            c = copy.deepcopy(self.case)
            c[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                v.validate(c)

    def test_dual_signs(self):
        for key, value in [('input_slope', '-1'), ('reset_slope', '1')]:
            c = copy.deepcopy(self.case)
            c[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                v.certify(c, enabled=True)

    def test_log_reciprocal_and_power_scaling(self):
        for q in [F(1), F(2), F(3, 5), F(10)**-107, F(10)**107]:
            a, b = v.log(q)
            c, d = v.log(1/q)
            self.assertLessEqual(a + c, 0)
            self.assertGreaterEqual(b + d, 0)
        with self.assertRaises(ValueError):
            v.log(F(0))

    def test_direct_fraction_series(self):
        def unit(q):
            u = (q-1)/(q+1)
            s, term = F(0), u
            for j in range(32):
                s += term/(2*j+1)
                term *= u*u
            return 2*s, 2*s + 2*term/(65*(1-u*u))
        l2, h2 = unit(F(2))
        p, t = v.validate(self.case)
        for q in [F(1), F(2), F(3, 5), F(2, 5), p[0], t[0], p[-1], t[-1]]:
            r, k = q, 0
            while r < 1:
                r *= 2
                k -= 1
            while r >= 2:
                r /= 2
                k += 1
            c, d = unit(r)
            c, d = (c+k*l2, d+k*h2) if k >= 0 else (c+k*h2, d+k*l2)
            a, b = v.log(q)
            self.assertLessEqual(max(a, c), min(b, d))

    def test_optimized_python_still_rejects(self):
        code = "import verify as v;c=v.load_case();c['p'][0]=str(int(c['p'][0])+1);v.validate(c)"
        r = subprocess.run([sys.executable, '-O', '-c', code], cwd=Path(__file__).parent, capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('ValueError', r.stderr)

    def test_cli_opt_in_and_no_writes(self):
        folder = Path(__file__).parent
        before = {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}
        r = subprocess.run([sys.executable, str(folder/'verify.py')], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        r = subprocess.run([sys.executable, str(folder/'verify.py'), '--enable'], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(json.loads(r.stdout)['compact_lower_decimal'], '4.9836437909')
        self.assertEqual(before, {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()})


if __name__ == '__main__':
    unittest.main()
