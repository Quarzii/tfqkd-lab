"""Local one-page laboratory interface for explicit TF-QKD inputs.

Run: python -m tfqkd.lab_web --port 8765
The server binds to loopback only and materializes uploaded CSV text in a
temporary directory before calling the existing laboratory wrapper.
"""
from __future__ import annotations

import argparse
import atexit
import html
import json
import math
import mimetypes
import secrets
import shutil
import tempfile
import threading
import time
import tomllib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlparse

from .lab_inputs import MissingInputs
from .measured_inputs import MeasurementRejected
from .web_inputs import FIELDS, RULE_MESSAGES, FieldValidationError, validate_web_input, get_path

MAX_REQUEST_BYTES = 4 * 1024 * 1024
MAX_SESSIONS = 8
SESSION_TTL_SECONDS = 60 * 60
ARTIFACT_ALLOW = {"report.md", "report.html", "result.json", "spectrum.png", "cumulative_variance.png"}
_SESSIONS: dict[str, dict] = {}
_LOCK = threading.Lock()
_CALCULATION_LOCK = threading.Lock()


def _cleanup_all():
    with _LOCK:
        sessions = list(_SESSIONS.values())
        _SESSIONS.clear()
    for entry in sessions:
        shutil.rmtree(entry["directory"], ignore_errors=True)


atexit.register(_cleanup_all)


def _prune_sessions(now=None):
    now = time.time() if now is None else now
    remove = []
    with _LOCK:
        for token, entry in _SESSIONS.items():
            if now - entry["created"] > SESSION_TTL_SECONDS:
                remove.append(token)
        if len(_SESSIONS) - len(remove) > MAX_SESSIONS:
            keep = max(0, MAX_SESSIONS - len(remove))
            active = [(entry["created"], token) for token, entry in _SESSIONS.items() if token not in remove]
            active.sort(reverse=True)
            retained = {token for _, token in active[:keep]}
            remove.extend(token for _, token in active if token not in retained)
        directories = []
        for token in remove:
            entry = _SESSIONS.pop(token, None)
            if entry:
                directories.append(entry["directory"])
    for directory in directories:
        shutil.rmtree(directory, ignore_errors=True)


def _asset(name: str) -> bytes:
    return (resources.files("tfqkd") / "web" / name).read_bytes()


def _safe_upload_name(name: str) -> str:
    candidate = PurePosixPath(str(name).replace("\\", "/")).name
    if not candidate or candidate in {".", ".."}:
        raise ValueError("Uploaded files need a plain filename")
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    if any(ch not in allowed for ch in candidate):
        raise ValueError(f"Unsafe upload filename: {name!r}")
    if not candidate.lower().endswith(".csv"):
        raise ValueError("Only CSV uploads are accepted")
    return candidate


def _safe_artifact_path(path_text: str) -> PurePosixPath:
    path = PurePosixPath(path_text)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError("Artifact path is not allowed")
    if len(path.parts) == 1 and path.name in ARTIFACT_ALLOW:
        return path
    if len(path.parts) == 2 and path.parts[0] in {"variant_1", "variant_2", "variant_3"} and path.parts[1] in ARTIFACT_ALLOW:
        return path
    raise ValueError("Artifact is not downloadable")


def _artifact_paths(directory: Path):
    paths = []
    for name in sorted(ARTIFACT_ALLOW):
        if (directory / name).exists():
            paths.append(name)
    for variant in ("variant_1", "variant_2", "variant_3"):
        for name in sorted(ARTIFACT_ALLOW):
            relative = PurePosixPath(variant) / name
            if (directory / relative).exists():
                paths.append(str(relative))
    return paths


def _walk_files(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "file":
                yield item
            yield from _walk_files(item)
    elif isinstance(value, list):
        for item in value:
            yield from _walk_files(item)


def _validate_toml_file_references(text: str, uploaded: set[str]):
    try:
        raw = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"TOML parse error: {error}") from error
    for reference in _walk_files(raw):
        if not isinstance(reference, str):
            raise ValueError("CSV file references must be strings")
        path = PurePosixPath(reference.replace("\\", "/"))
        if path.is_absolute() or ".." in path.parts or len(path.parts) != 1:
            raise ValueError(f"CSV file reference is not a local upload name: {reference}")
        if path.name not in uploaded:
            raise ValueError(f"CSV file was referenced but not uploaded: {path.name}")
    return raw


def _materialize_files(directory: Path, files):
    uploaded = set()
    if files is not None and not isinstance(files, list):
        raise ValueError("Uploaded files must be a list")
    for item in files or []:
        if not isinstance(item, dict):
            raise ValueError("Each uploaded file must be an object")
        if set(item) - {"name", "content"}:
            raise ValueError("Uploaded file objects accept only name and content")
        name = _safe_upload_name(item.get("name", ""))
        if name in uploaded:
            raise ValueError(f"Duplicate upload filename: {name}")
        content = item.get("content", "")
        if not isinstance(content, str):
            raise ValueError(f"Upload {name} content must be text")
        (directory / name).write_text(content)
        uploaded.add(name)
    return uploaded


def _canonical_toml(raw):
    """Serialize parsed web configuration; only representation changes."""
    def scalar(item):
        if isinstance(item, dict):
            return '{ '+', '.join(json.dumps(k)+' = '+scalar(v) for k,v in item.items())+' }'
        if isinstance(item, list):return '['+', '.join(scalar(v) for v in item)+']'
        if isinstance(item,(str,int,float,bool)):return json.dumps(item,ensure_ascii=False)
        raise ValueError('Configuration values must be numbers, strings, booleans or arrays/tables')
    if any(not isinstance(values,dict) for values in raw.values()):raise ValueError("Use named TOML sections for configuration inputs")
    return '\n'.join('['+json.dumps(section)+']\n'+'\n'.join(json.dumps(k)+' = '+scalar(v) for k,v in values.items())
                     for section,values in raw.items())+'\n'


def _web_context(raw, payload):
    warnings={}
    line=raw.get('line',{})
    ignored=[]
    if 'loss_a_db' in line and 'loss_b_db' in line:
        for section in ('line','keyrate'):
            if 'attenuation_db_per_km' in raw.get(section,{}):
                ignored.append(section+'.attenuation_db_per_km')
                raw[section].pop('attenuation_db_per_km')
        if ignored:warnings['line.attenuation_db_per_km']='Complete arm losses take precedence; per-km attenuation was not used.'
    warnings.update(validate_web_input(raw))
    intensity_rows=[]
    p=raw.get('keyrate',{})
    # Appendix D, bertaina2024 / QKD.ipynb Cell 21: SNS decoys use two-user intensities.
    # Appendix B, bertaina2024 / QKD.ipynb Cell 19: CAL alpha^2 is already per user.
    keys=('decoy_big','decoy_medium','decoy_mini') if raw.get('protocol',{}).get('name')=='SNS-AOPP' else ('u_cal',)
    supplied=payload.get('intensities_per_user')
    if supplied is not None and not isinstance(supplied,dict):raise ValueError('intensities_per_user must be an object')
    for key in keys:
        if key not in p:continue
        factor=2 if key.startswith('decoy_') else 1
        per_user=p[key]/factor
        if supplied is not None and key in supplied:
            try:matches=math.isclose(float(supplied[key])*factor,p[key],rel_tol=1e-12,abs_tol=0.)
            except (ValueError,TypeError):matches=False
            if not matches:raise FieldValidationError({'keyrate.'+key:'Per-user input does not match the core configuration; rebuild it or remove form metadata for an Advanced edit.'})
        intensity_rows.append(dict(parameter=key,per_user=per_user,core=p[key],factor=factor,
            source='Appendix D / QKD.ipynb Cell 21, bertaina2024' if factor==2 else 'Appendix B / QKD.ipynb Cell 19, bertaina2024'))
    return dict(warnings=warnings,ignored_attenuation=ignored,intensities=intensity_rows,
                source=payload.get('example_source'),entry='structured form' if supplied is not None else 'Advanced core TOML')


def _example(name):
    from .config import ROOT
    paths={'zhou2023':ROOT/'examples/measured_inputs/zhou2023_403.73_ideal.toml',
           'bertaina2024':ROOT/'examples/bertaina2024_table3.toml'}
    if name not in paths:raise ValueError('Unknown published example')
    raw=tomllib.loads(paths[name].read_text());values={}
    for id,spec in FIELDS.items():
        if not spec['path']:continue
        v=get_path(raw,spec['path'])
        key=spec['path'].split('.')[-1]
        if v is None and spec['path'].startswith('laser.'):v=raw.get('physics',{}).get(key)
        if v is None and id in ('fiberL','fc1'):v=raw.get('physics',{}).get(key)
        if v is None and id in ('alpha','detectorEfficiency','detectorDark','detectorError'):
            alias={'alpha':'attenuation_db_per_km','detectorEfficiency':'detector_efficiency','detectorDark':'detector_dark_count_rate_hz','detectorError':'detector_error'}[id]
            v=raw.get('keyrate',{}).get(alias)
        if v is not None:
            # Appendix D, bertaina2024: displayed SNS values are per user.
            values[id]=v/2 if id in ('decoyBig','decoyMedium','decoyMini') else v
    if name=='zhou2023':
        values.update(mode='phase',protocol='SNS-AOPP',detectorMode='two',reach='false',length=403.73,imbalance=.01)
        source='Zhou2023 Tables S1–S5, Methods (f_EC), Fig.3e–g; per-user SNS intensities 0.493 / 0.105 / 0.0002. Quantum-channel phase RMS and intrinsic e_d at 403 km are not published: sigma=0 is an explicit ideal upper-phase bound, not a measurement. tau=200 ns and tau_PS=0 retain the approved timing idealization; the 500 MHz effective clock already excludes reference slots.'
    else:
        values.update(mode='spectral',protocol=raw['protocol']['name'],detectorMode='scalar',reach='false',laserInput='coefficients',actuatorInput='omega',laserModel=raw['laser']['model'],schemeLasers=raw['scheme']['lasers'],schemeCompensation=raw['scheme']['compensation'])
        source='Bertaina2024 Table III (noise), Table II / Appendix D and QKD.ipynb (protocol), Table I demonstration geometry. The actuator pole is an explicit engineering example, not a measured published specification.'
    return dict(ok=True,fields=values,source=source,configuration=paths[name].read_text())


def _write_single_input(payload, directory: Path) -> Path:
    text = payload.get("toml")
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Supply TOML text for the calculation")
    uploaded = _materialize_files(directory, payload.get("files"))
    raw = _validate_toml_file_references(text, uploaded)
    context = _web_context(raw, payload)
    path = directory / "input.toml"
    path.write_text(_canonical_toml(raw))
    (directory / "web_context.json").write_text(json.dumps(context))
    return path


def _write_compare_input(payload, directory: Path) -> Path:
    variants = payload.get("variants")
    if not isinstance(variants, list) or len(variants) not in (2, 3):
        raise ValueError("Comparison needs two or three variants")
    labels, paths = [], []
    for index, variant in enumerate(variants, start=1):
        if not isinstance(variant, dict):
            raise ValueError("Each comparison variant must be an object")
        label = str(variant.get("label", "")).strip()
        if not label:
            raise ValueError("Each comparison variant needs a label")
        label = label[:80]
        variant_dir = directory / f"variant_{index}_input"
        variant_dir.mkdir()
        try:path = _write_single_input(variant, variant_dir)
        except FieldValidationError as error:
            raise FieldValidationError({f'variants.{index-1}.'+key:message for key,message in error.field_errors.items()}) from error
        labels.append(label)
        paths.append(str(path))
    if len(set(labels)) != len(labels):
        raise ValueError("Comparison labels must be distinct")
    lengths = payload.get("working_total_lengths_km")
    if not isinstance(lengths, list) or not lengths:
        raise ValueError("working_total_lengths_km must be a nonempty list")
    if any(not isinstance(x,(int,float)) or isinstance(x,bool) or not math.isfinite(x) or x<=0 for x in lengths):
        raise FieldValidationError({'comparison.working_total_lengths_km':'All working total lengths must be finite and > 0 [km].'})
    workers = payload.get("workers", 1)
    if type(workers) is not int or workers != 1:
        raise ValueError("The local web server runs comparison workers=1")
    body = ["[comparison]"]
    body.append("configurations = [" + ", ".join(json.dumps(path) for path in paths) + "]")
    body.append("labels = [" + ", ".join(json.dumps(label) for label in labels) + "]")
    body.append("working_total_lengths_km = [" + ", ".join(str(float(x)) for x in lengths) + "]")
    body.append("workers = 1")
    path = directory / "comparison.toml"
    path.write_text("\n".join(body) + "\n")
    return path


def _session(directory: Path, markdown: str, result: dict, kind: str):
    token = secrets.token_urlsafe(18)
    with _LOCK:
        _SESSIONS[token] = dict(directory=directory, created=time.time(), kind=kind)
    _prune_sessions()
    return dict(
        ok=True,
        kind=kind,
        markdown=markdown,
        result=result,
        artifacts={
            name: f"/artifacts/{token}/{name}"
            for name in _artifact_paths(directory)
        },
    )


def _run_payload(payload):
    directory = Path(tempfile.mkdtemp(prefix="tfqkd-web-"))
    try:
        config = _write_single_input(payload, directory)
        from .lab_report import write_report
        from .lab_run import markdown, run

        with _CALCULATION_LOCK:
            result = run(config, workers=1)
            result["web_input"] = json.loads((directory / "web_context.json").read_text())
            write_report(result, directory)
        return _session(directory, markdown(result), result, "single")
    except Exception:
        shutil.rmtree(directory, ignore_errors=True)
        raise


def _compare_payload(payload):
    directory = Path(tempfile.mkdtemp(prefix="tfqkd-web-"))
    try:
        config = _write_compare_input(payload, directory)
        from .lab_compare import compare, comparison_markdown, write_comparison

        with _CALCULATION_LOCK:
            result = compare(config, workers=1)
            for index, variant in enumerate(result['variants'],start=1):
                context = directory / f'variant_{index}_input' / 'web_context.json'
                variant['result']['web_input'] = json.loads(context.read_text())
            write_comparison(result, directory)
        return _session(directory, comparison_markdown(result), result, "comparison")
    except Exception:
        shutil.rmtree(directory, ignore_errors=True)
        raise


def _json_response(handler, status: int, payload: dict):
    data = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(data)))
    handler.send_header("Cache-Control", "no-store")
    handler.end_headers()
    handler.wfile.write(data)


def _error_payload(error: Exception):
    return dict(
        ok=False,
        error=html.escape(str(error)),
        missing_inputs=getattr(error, "parameters", None),
        diagnostics=getattr(error, "diagnostics", None),
        field_errors=getattr(error, "field_errors", None),
    )


class LabWebHandler(BaseHTTPRequestHandler):
    server_version = "TFQKDLabWeb/1.0"

    def log_message(self, fmt, *args):
        return

    def _authorized(self):
        host = self.headers.get("Host", "")
        allowed_hosts = {
            f"127.0.0.1:{self.server.server_port}",
            f"localhost:{self.server.server_port}",
            f"[::1]:{self.server.server_port}",
        }
        if host not in allowed_hosts:
            return False
        origin = self.headers.get("Origin")
        if origin:
            try:
                parsed = urlparse(origin)
                origin_host = f"{parsed.hostname}:{parsed.port}"
            except ValueError:
                return False
            if origin_host not in {h.replace("[", "").replace("]", "") for h in allowed_hosts}:
                return False
        client = self.client_address[0]
        return client in {"127.0.0.1", "::1"} or client.startswith("::ffff:127.")

    def _send_asset(self, name: str, content_type: str):
        data = _asset(name)
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        _prune_sessions()
        parsed = urlparse(self.path)
        if not self._authorized():
            _json_response(self, 403, dict(ok=False, error="Local loopback Host/Origin required"))
            return
        if parsed.path == "/":
            self._send_asset("index.html", "text/html; charset=utf-8")
            return
        if parsed.path == "/styles.css":
            self._send_asset("styles.css", "text/css; charset=utf-8")
            return
        if parsed.path == "/app.js":
            self._send_asset("app.js", "text/javascript; charset=utf-8")
            return
        if parsed.path == "/api/schema":
            _json_response(self, 200, dict(ok=True, fields=FIELDS, rules=RULE_MESSAGES))
            return
        if parsed.path.startswith("/api/examples/"):
            try:_json_response(self, 200, _example(parsed.path.rsplit('/',1)[-1]))
            except ValueError as error:_json_response(self, 404, _error_payload(error))
            return
        if parsed.path.startswith("/artifacts/"):
            parts = [unquote(part) for part in parsed.path.split("/") if part]
            if len(parts) < 3:
                _json_response(self, 404, dict(ok=False, error="Artifact path not found"))
                return
            _, token, *artifact_parts = parts
            try:
                relative = _safe_artifact_path("/".join(artifact_parts))
            except ValueError as error:
                _json_response(self, 404, dict(ok=False, error=str(error)))
                return
            with _LOCK:
                entry = _SESSIONS.get(token)
            if not entry:
                _json_response(self, 404, dict(ok=False, error="Artifact session expired"))
                return
            root = Path(entry["directory"]).resolve()
            path = (root / relative).resolve()
            if not path.is_relative_to(root) or not path.exists():
                _json_response(self, 404, dict(ok=False, error="Artifact is missing"))
                return
            data = path.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", mimetypes.guess_type(path.name)[0] or "application/octet-stream")
            self.send_header("Content-Length", str(len(data)))
            disposition = "inline" if path.suffix == ".html" else "attachment"
            self.send_header("Content-Disposition", f'{disposition}; filename="{path.name}"')
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)
            return
        _json_response(self, 404, dict(ok=False, error="Not found"))

    def do_POST(self):
        _prune_sessions()
        if not self._authorized():
            _json_response(self, 403, dict(ok=False, error="Local loopback Host/Origin required"))
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            _json_response(self, 400, dict(ok=False, error="Invalid Content-Length"))
            return
        if length <= 0 or length > MAX_REQUEST_BYTES:
            _json_response(self, 413, dict(ok=False, error="Request body is empty or exceeds 4 MiB"))
            return
        try:
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("JSON request body must be an object")
            if self.path == "/api/run":
                response = _run_payload(payload)
            elif self.path == "/api/compare":
                response = _compare_payload(payload)
            else:
                _json_response(self, 404, dict(ok=False, error="Not found"))
                return
            _json_response(self, 200, response)
        except (MissingInputs, MeasurementRejected, ValueError, TypeError, ArithmeticError, OSError) as error:
            _json_response(self, 422, _error_payload(error))


def serve(port: int, host: str = "127.0.0.1"):
    if host not in {"127.0.0.1", "::1", "localhost"}:
        raise ValueError("The laboratory web interface only binds to loopback")
    server = ThreadingHTTPServer((host, port), LabWebHandler)
    actual = server.server_address
    print(f"TF-QKD lab web interface: http://{actual[0]}:{actual[1]}/")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        _cleanup_all()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()
    serve(args.port, args.host)


if __name__ == "__main__":
    main()
