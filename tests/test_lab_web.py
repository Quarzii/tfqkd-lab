import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest import mock

from tfqkd.config import ROOT
import tfqkd.lab_web as lab_web
from tfqkd.lab_web import LabWebHandler


class NonClosingBytesIO(BytesIO):
    def close(self):
        pass


class FakeSocket:
    def __init__(self, request):
        self.input = BytesIO(request)
        self.output = NonClosingBytesIO()

    def makefile(self, mode, buffering=None):
        return self.input if "r" in mode else self.output

    def sendall(self, data):
        self.output.write(data)


class FakeServer:
    server_port = 8123


def request(method, path, body=None, host=None, extra_headers=None):
    headers = {"Host": host or "127.0.0.1:8123"}
    headers.update(extra_headers or {})
    payload = None
    if body is not None:
        payload = json.dumps(body).encode()
        headers["Content-Type"] = "application/json"
        headers["Content-Length"] = str(len(payload))
    lines = [f"{method} {path} HTTP/1.1"]
    lines += [f"{key}: {value}" for key, value in headers.items()]
    raw = ("\r\n".join(lines) + "\r\n\r\n").encode() + (payload or b"")
    sock = FakeSocket(raw)
    LabWebHandler(sock, ("127.0.0.1", 43210), FakeServer())
    data = sock.output.getvalue()
    head, _, body_bytes = data.partition(b"\r\n\r\n")
    header_lines = head.decode("iso-8859-1").split("\r\n")
    status = int(header_lines[0].split()[1])
    parsed = []
    for line in header_lines[1:]:
        if ":" in line:
            key, value = line.split(":", 1)
            parsed.append((key, value.strip()))
    return status, parsed, body_bytes


class LabWebTests(unittest.TestCase):
    def test_measured_phase_run_writes_downloadable_report(self):
        toml = (ROOT / "examples/measured_inputs/zhou2023_403.73_ideal.toml").read_text()
        status, _, data = request("POST", "/api/run", {"toml": toml, "files": []})
        self.assertEqual(status, 200)
        payload = json.loads(data)
        self.assertEqual(payload["kind"], "single")
        self.assertIn("report.html", payload["artifacts"])
        self.assertIn("residual phase", payload["markdown"].lower())
        html_status, headers, html_data = request("GET", payload["artifacts"]["report.html"])
        self.assertEqual(html_status, 200)
        self.assertIn(b"Twin-Field QKD laboratory report", html_data)
        disposition = next(v for k, v in headers if k.lower() == "content-disposition")
        self.assertTrue(disposition.startswith("inline"))
        json_status, json_headers, _ = request("GET", payload["artifacts"]["result.json"])
        self.assertEqual(json_status, 200)
        json_disposition = next(v for k, v in json_headers if k.lower() == "content-disposition")
        self.assertTrue(json_disposition.startswith("attachment"))

    def test_rejects_nonlocal_csv_reference_before_core_read(self):
        toml = """
[laser]
model = "free"
spectrum = { file = "../secret.csv", mode = "direct", quantity = "phase", frequency_unit = "Hz", psd_unit = "rad^2/Hz", sidedness = "one-sided", pass = "single" }
"""
        status, _, data = request("POST", "/api/run", {"toml": toml, "files": []})
        self.assertEqual(status, 422)
        payload = json.loads(data)
        self.assertIn("not a local upload name", payload["error"])

    def test_ui_has_no_apparatus_or_phase_defaults(self):
        index = (ROOT / "tfqkd/web/index.html").read_text()
        for forbidden in ('value="0.05"', 'value="2e-7"', 'value="0"', "efficiency = 0.6", "clockrate_hz = 500000000.0"):
            self.assertNotIn(forbidden, index)
        self.assertIn("Estimate key rate for your optical link.", index)
        script = (ROOT / "tfqkd/web/app.js").read_text()
        self.assertNotIn('name = ""', script)
        for control in ("laserModel", "laserInput", "laserCsvMode", "lineCsvMode", "schemeCompensation",
                        "actuatorInput", "actuatorCsvFile", "sigmaLimit", "clockrate", "detectorEfficiency"):
            self.assertIn(f'id="{control}"', index)
        self.assertIn("JSON.stringify", script)

    def test_escapes_error_text_and_requires_loopback_host(self):
        bad_host_status, _, _ = request("GET", "/", host="example.com")
        self.assertEqual(bad_host_status, 403)
        bad_origin_status, _, _ = request("GET", "/", extra_headers={"Origin": "http://127.0.0.1:notaport"})
        self.assertEqual(bad_origin_status, 403)
        toml = """
[laser]
model = "free"
spectrum = { file = "evil<script>.csv", mode = "direct", quantity = "phase", frequency_unit = "Hz", psd_unit = "rad^2/Hz", sidedness = "one-sided", pass = "single" }
"""
        status, _, data = request("POST", "/api/run", {"toml": toml, "files": []})
        self.assertEqual(status, 422)
        payload = json.loads(data)
        self.assertNotIn("<script>", payload["error"])
        self.assertIn("&lt;script&gt;", payload["error"])
        bad_files_status, _, bad_files_data = request("POST", "/api/run", {"toml": "[phase]\nmode='measured'\n", "files": "x"})
        self.assertEqual(bad_files_status, 422)
        self.assertIn("Uploaded files must be a list", json.loads(bad_files_data)["error"])

    def test_sessions_are_pruned_to_finite_retention(self):
        lab_web._cleanup_all()
        directories = [Path(tempfile.mkdtemp(prefix="tfqkd-web-test-")) for _ in range(lab_web.MAX_SESSIONS + 2)]
        try:
            for directory in directories:
                (directory / "report.md").write_text("# x\n")
                lab_web._session(directory, "# x\n", {"ok": True}, "single")
            self.assertLessEqual(len(lab_web._SESSIONS), lab_web.MAX_SESSIONS)
            self.assertFalse(directories[0].exists())
        finally:
            lab_web._cleanup_all()
            for directory in directories:
                if directory.exists():
                    lab_web.shutil.rmtree(directory, ignore_errors=True)

    def test_compare_endpoint_materializes_variants_and_artifacts(self):
        result = {"mode": "comparison", "variants": [], "working_lengths": [], "input_case_calculations": 0, "spectral_calculations": 0}

        def fake_write(_, directory):
            Path(directory).mkdir(parents=True, exist_ok=True)
            (Path(directory) / "report.md").write_text("# comparison\n")
            (Path(directory) / "report.html").write_text("<!doctype html><p>comparison</p>")
            (Path(directory) / "result.json").write_text(json.dumps(result))
            variant = Path(directory) / "variant_1"
            variant.mkdir()
            (variant / "report.html").write_text("<!doctype html><p>variant</p>")
            (variant / "report.md").write_text("# variant\n")
            (variant / "result.json").write_text(json.dumps({"variant": 1}))

        with mock.patch("tfqkd.lab_compare.compare", return_value=result) as compare_mock, \
             mock.patch("tfqkd.lab_compare.comparison_markdown", return_value="# comparison\n"), \
             mock.patch("tfqkd.lab_compare.write_comparison", side_effect=fake_write):
            body = {
                "variants": [
                    {"label": "A", "toml": "[phase]\nmode='measured'\nsigma_phi_rad=0\ntau_s=1e-6\ntau_ps_s=0\n", "files": []},
                    {"label": "B", "toml": "[phase]\nmode='measured'\nsigma_phi_rad=0.1\ntau_s=1e-6\ntau_ps_s=0\n", "files": []},
                ],
                "working_total_lengths_km": [50.0, 100.0],
                "workers": 1,
            }
            status, _, data = request("POST", "/api/compare", body)
            self.assertEqual(status, 200)
            payload = json.loads(data)
            self.assertEqual(payload["kind"], "comparison")
            self.assertIn("result.json", payload["artifacts"])
            self.assertIn("variant_1/report.html", payload["artifacts"])
            html_status, headers, html_data = request("GET", payload["artifacts"]["variant_1/report.html"])
            self.assertEqual(html_status, 200)
            self.assertIn(b"variant", html_data)
            disposition = next(v for k, v in headers if k.lower() == "content-disposition")
            self.assertTrue(disposition.startswith("inline"))
            token = payload["artifacts"]["variant_1/report.html"].split("/")[2]
            traversal_status, _, _ = request("GET", f"/artifacts/{token}/../report.html")
            self.assertEqual(traversal_status, 404)
            compare_mock.assert_called_once()
            body["workers"] = None
            bad_status, _, bad_data = request("POST", "/api/compare", body)
            self.assertEqual(bad_status, 422)
            self.assertIn("workers=1", json.loads(bad_data)["error"])


if __name__ == "__main__":
    unittest.main()
