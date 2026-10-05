const $ = (id) => document.getElementById(id);
let fields = {}, ruleMessages = {}, editorDirty = false, exampleSource = "";
let builtText = "";

function enabledField(id) {
  const el = $(id);
  if (!el || el.disabled || ["labelA","labelB","labelC","tomlA","tomlB","tomlC","lengths","toml"].includes(id)) return false;
  if (el.closest(".spectral-only") && value("mode") !== "spectral") return false;
  if (el.closest(".phase-fields") && value("mode") !== "phase") return false;
  if (el.closest(".protocol-cal") && value("protocol") !== "CAL") return false;
  if (el.closest(".protocol-sns") && value("protocol") === "CAL") return false;
  if (el.closest(".detector-two") && value("detectorMode") !== "two") return false;
  if (el.closest(".detector-scalar") && value("detectorMode") !== "scalar") return false;
  if (el.closest(".cavity-fields") && value("laserModel") !== "cavity") return false;
  if (el.closest(".dual-fields") && value("schemeCompensation") !== "dual") return false;
  return true;
}

function fieldMessage(id, text, severity = "error") {
  const el = $(id), message = $(`${id}Message`);
  if (!el || !message) return;
  message.textContent = text;
  if (text && severity === "error" && el.closest("details")) el.closest("details").open = true;
  message.className = `field-message field-${severity}`;
  el.setAttribute("aria-invalid", text && severity === "error" ? "true" : "false");
}

function clearMessages() {
  for (const id of Object.keys(fields)) fieldMessage(id, "");
}

function rulePass(n, rule) {
  if (!Number.isFinite(n)) return false;
  return { positive: () => n > 0, nonnegative: () => n >= 0, probability: () => n > 0 && n <= 1, sns_probability: () => n > 0 && n < 1,
    error: () => n >= 0 && n < 0.5, phase: () => n >= 0 && n <= Math.PI,
    phase_positive: () => n > 0 && n <= Math.PI, fec: () => n >= 1,
    integer: () => n >= 0 && Number.isInteger(n), correlation: () => n === 2 || n === 4,
    finite: () => true }[rule]();
}

function validateForm() {
  clearMessages();
  let valid = true;
  for (const [id, spec] of Object.entries(fields)) {
    if (!spec.rule || !enabledField(id)) continue;
    const el = $(id), text = value(id);
    if (!text && !el.validity.badInput) continue;
    const n = Number(text);
    if (el.validity.badInput || !rulePass(n, spec.rule)) {
      fieldMessage(id, ruleMessages[spec.rule]); valid = false; continue;
    }
    let warning = "";
    if (spec.unusual === "high_phase" && n > 1) warning = "Unusually large phase RMS (> 1 rad); calculation is allowed.";
    if (spec.unusual === "low_efficiency" && n < 0.05) warning = "Unusually low detector efficiency (< 0.05); calculation is allowed.";
    if (spec.unusual === "high_fec" && n > 2) warning = "Unusually large f_EC (> 2); calculation is allowed.";
    if (spec.unusual === "sns_typical" && (n < 0.05 || n > 0.3)) warning = "Outside the typical sending range 0.05-0.3; calculation is allowed.";
    fieldMessage(id, warning, "warning");
  }
  if (value("protocol") !== "CAL" && ["decoyBig","decoyMedium","decoyMini"].every((id) => value(id))) {
    const [big, medium, mini] = ["decoyBig","decoyMedium","decoyMini"].map((id) => Number(value(id)));
    if (!(big > medium && medium > mini && mini >= 0)) {
      for (const id of ["decoyBig","decoyMedium","decoyMini"]) fieldMessage(id, "Require big > medium > mini >= 0.");
      valid = false;
    }
  }
  for (const [id, other] of [["lossA","lossB"],["lossB","lossA"]]) {
    if (value(id) && !value(other)) { fieldMessage(other, "Supply both complete arm losses."); valid = false; }
  }
  if (value("length") && value("imbalance") && Math.abs(Number(value("imbalance"))) >= Number(value("length"))) {
    fieldMessage("imbalance", "Magnitude must be less than total length (both arms positive)."); valid = false;
  }
  return valid;
}

function updateLossPriority() {
  const complete = value("lossA") !== "" && value("lossB") !== "";
  $("alpha").disabled = complete;
  if (complete) fieldMessage("alpha", "Complete arm losses take precedence; attenuation per km is not used.", "warning");
}

function showFieldErrors(errors) {
  for (const [path, text] of Object.entries(errors || {})) {
    const variant = path.match(/^variants\.([012])\.(.*)$/);
    if (variant) { fieldMessage("toml" + "ABC"[Number(variant[1])], variant[2] + ": " + text); continue; }
    if (path === "comparison.working_total_lengths_km") { fieldMessage("lengths", text); continue; }
    for (const [id, spec] of Object.entries(fields)) {
      let match = spec.path === path;
      const aliases = { "keyrate.detector_error": value("detectorMode") === "two" ? "twoDetectorError" : "detectorError",
        "keyrate.detector_efficiency": "detectorEfficiency", "keyrate.detector_dark_count_rate_hz": "detectorDark",
        "keyrate.attenuation_db_per_km": "alpha" };
      if (aliases[path] === id) match = true;
      if (path.startsWith("physics.") && ["laser.","line."].some((prefix) => spec.path === prefix + path.split(".").at(-1))) match = true;
      if (path === "protocol.name" && id === "protocol") match = true;
      if (path === "scheme.lasers" && id === "schemeLasers") match = true;
      if (path === "scheme.compensation" && id === "schemeCompensation") match = true;
      if (match) fieldMessage(id, text);
    }
  }
}

function decorateFields() {
  for (const [id, spec] of Object.entries(fields)) {
    const el = $(id); if (!el) continue;
    const label = el.closest("label"); if (!label) continue;
    for (const child of Array.from(label.childNodes)) if (child.nodeType === Node.TEXT_NODE) child.remove();
    const heading = document.createElement("span"); heading.className = "field-heading";
    heading.append(document.createTextNode(spec.label));
    const tip = document.createElement("button"); tip.type = "button"; tip.className = "field-help"; tip.textContent = "?";
    tip.setAttribute("aria-label", `${spec.label}: ${spec.help}`);
    function showTip() {
      const box = $("fieldTooltip"); box.textContent = spec.help; box.hidden = false;
      const r = tip.getBoundingClientRect();
      box.style.left = `${Math.max(12, Math.min(r.left, innerWidth - box.offsetWidth - 12))}px`;
      box.style.top = `${Math.max(12, r.bottom + box.offsetHeight + 8 > innerHeight ? r.top - box.offsetHeight - 8 : r.bottom + 8)}px`;
    }
    tip.addEventListener("mouseenter", showTip); tip.addEventListener("focus", showTip);
    tip.addEventListener("mouseleave", () => $("fieldTooltip").hidden = true);
    tip.addEventListener("blur", () => $("fieldTooltip").hidden = true);
    heading.append(tip);
    const variable = document.createElement("small"); variable.className = "field-variable"; variable.textContent = spec.variable;
    const message = document.createElement("span"); message.id = `${id}Message`; message.className = "field-message"; message.setAttribute("aria-live", "polite");
    label.prepend(heading, variable); label.append(message); el.setAttribute("aria-describedby", message.id);
  }
}

async function loadExample() {
  const response = await fetch(`/api/examples/${value("exampleChoice")}`);
  const data = await response.json(); if (!response.ok) throw new Error(data.error);
  for (const id of Object.keys(fields)) {
    const el = $(id); if (!el || ["exampleChoice","files","labelA","labelB","labelC","tomlA","tomlB","tomlC","lengths"].includes(id)) continue;
    el.value = "";
  }
  for (const [id, val] of Object.entries(data.fields)) if ($(id)) $(id).value = String(val);
  exampleSource = data.source; $("exampleSource").textContent = data.source;
  updateMode(); updateLossPriority(); buildToml();
  showStatus("Example loaded. Source values and explicit phase/timing assumptions are shown above; review before running.");
}

function value(id) {
  return $(id).value.trim();
}

function tomlString(text) {
  return JSON.stringify(text);
}

function finiteNumber(id, label) {
  const text = value(id);
  if (!text) return null;
  const number = Number(text);
  if (!Number.isFinite(number)) throw new Error(`${label} must be a finite number`);
  return text;
}

function numberLine(key, id, label) {
  const text = finiteNumber(id, label || key);
  return text ? `${key} = ${text}` : "";
}

function textLine(key, id) {
  const text = value(id);
  return text ? `${key} = ${tomlString(text)}` : "";
}

function section(name, lines) {
  const body = lines.filter(Boolean).join("\n");
  return body ? `[${name}]\n${body}\n` : "";
}

function tableArray(entries) {
  const pairs = entries.filter(([, v]) => v !== null && v !== "" && v !== undefined);
  if (!pairs.length) return "";
  return "{ " + pairs.map(([k, v]) => `${k} = ${v}`).join(", ") + " }";
}

function quotedValue(id) {
  const text = value(id);
  return text ? tomlString(text) : null;
}

function numberValue(id, label) {
  const text = finiteNumber(id, label);
  return text || null;
}

function updateMode() {
  const mode = value("mode");
  const protocol = value("protocol");
  const detector = value("detectorMode");
  const compensation = value("schemeCompensation");
  document.body.classList.toggle("spectral", mode === "spectral");
  document.body.classList.toggle("protocol-cal-selected", protocol === "CAL");
  document.body.classList.toggle("detector-two-selected", detector === "two");
  document.body.classList.toggle("dual-selected", compensation === "dual");
  document.body.classList.toggle("cavity-selected", value("laserModel") === "cavity");
}

function csvSpec(prefix, node) {
  const file = value(`${prefix}CsvFile`);
  if (!file) return "";
  const entries = [
    ["file", tomlString(file)],
    ["mode", quotedValue(`${prefix}CsvMode`)],
    ["quantity", quotedValue(`${prefix}CsvQuantity`)],
    ["frequency_unit", tomlString("Hz")],
    ["psd_unit", quotedValue(`${prefix}CsvPsdUnit`)],
    ["sidedness", quotedValue(`${prefix}CsvSidedness`)],
    ["pass", quotedValue(`${prefix}CsvPass`)],
    ["round_trip_psd_factor", numberValue(`${prefix}RoundTripFactor`, `${node} round-trip PSD factor`)],
    ["extrapolation", quotedValue(`${prefix}Extrapolation`)]
  ];
  if (node === "line") entries.push(["measurement_length_km", numberValue("lineCsvLength", "line measurement length")]);
  const result = tableArray(entries);
  if (!result.includes("mode =") || !result.includes("quantity =") || !result.includes("psd_unit =") || !result.includes("sidedness =") || !result.includes("pass =")) {
    throw new Error(`${node} CSV needs mode, quantity, PSD unit, sidedness and pass`);
  }
  if (node === "line" && !result.includes("measurement_length_km =")) {
    throw new Error("line CSV needs measured length");
  }
  return result;
}

function actuatorResponseSpec() {
  const file = value("actuatorCsvFile");
  if (!file) return "";
  const result = tableArray([
    ["file", tomlString(file)],
    ["frequency_unit", tomlString("Hz")],
    ["magnitude_unit", quotedValue("actuatorMagnitudeUnit")],
    ["phase_unit", quotedValue("actuatorPhaseUnit")]
  ]);
  if (!result.includes("magnitude_unit =") || !result.includes("phase_unit =")) {
    throw new Error("actuator response CSV needs magnitude and phase units");
  }
  return result;
}

function detectorSection() {
  if (value("detectorMode") === "two") {
    const d0 = tableArray([
      ["name", tomlString("D0")],
      ["efficiency", numberValue("d0Efficiency", "D0 efficiency")],
      ["dark_count_rate_hz", numberValue("d0Dark", "D0 dark count")],
      ["background_count_rate_hz", numberValue("d0Background", "D0 background")]
    ]);
    const d1 = tableArray([
      ["name", tomlString("D1")],
      ["efficiency", numberValue("d1Efficiency", "D1 efficiency")],
      ["dark_count_rate_hz", numberValue("d1Dark", "D1 dark count")],
      ["background_count_rate_hz", numberValue("d1Background", "D1 background")]
    ]);
    const lines = [];
    if (d0 || d1) lines.push(`channels = [${d0 || "{ name = \"D0\" }"}, ${d1 || "{ name = \"D1\" }"}]`);
    lines.push(numberLine("error", "twoDetectorError", "detector error"));
    return section("detector", lines);
  }
  return section("detector", [
    numberLine("efficiency", "detectorEfficiency", "detector efficiency"),
    numberLine("dark_count_rate_hz", "detectorDark", "detector dark count"),
    numberLine("error", "detectorError", "detector error")
  ]);
}

function intensityLine(key, id) {
  const text = finiteNumber(id, key);
  // Appendix D, bertaina2024 / QKD.ipynb Cell 21: total SNS intensity is twice the per-user input.
  return text ? `${key} = ${2 * Number(text)}` : "";
}

function buildToml() {
  updateLossPriority();
  if (!validateForm()) { updateLossPriority(); throw new Error("Correct the highlighted fields before building or calculating."); }
  updateLossPriority();
  const mode = value("mode");
  const parts = [];
  if (mode === "phase") {
    parts.push(section("phase", [
      'mode = "measured"',
      numberLine("sigma_phi_rad", "sigma", "sigma_phi"),
      numberLine("tau_s", "tau", "tau"),
      numberLine("tau_ps_s", "taups", "tau_PS")
    ]));
  }
  const lineLines = [
    numberLine("length_km", "length", "total length"),
    numberLine("imbalance_km", "imbalance", "imbalance"),
    numberLine("loss_a_db", "lossA", "loss A"),
    numberLine("loss_b_db", "lossB", "loss B"),
    $("alpha").disabled ? "" : numberLine("attenuation_db_per_km", "alpha", "attenuation"),
    mode === "spectral" ? numberLine("l", "fiberL", "fiber l") : "",
    mode === "spectral" ? numberLine("fc1_hz", "fc1", "fiber fc1") : ""
  ];
  parts.push(detectorSection());
  parts.push(section("keyrate", [
    numberLine("clockrate_hz", "clockrate", "clock"),
    numberLine("f_error", "fError", "f_EC"),
    value("protocol") !== "CAL" ? intensityLine("decoy_big", "decoyBig") : "",
    value("protocol") !== "CAL" ? intensityLine("decoy_medium", "decoyMedium") : "",
    value("protocol") !== "CAL" ? intensityLine("decoy_mini", "decoyMini") : "",
    value("protocol") !== "CAL" ? numberLine("pz_sns", "pzSns", "SNS p_z") : "",
    value("protocol") !== "CAL" ? numberLine("eps_sns_aopp", "epsSns", "SNS epsilon") : "",
    value("protocol") === "CAL" ? numberLine("pz_cal", "pzCal", "CAL p_z") : "",
    value("protocol") === "CAL" ? numberLine("nmin_cal", "nminCal", "CAL n_min") : "",
    value("protocol") === "CAL" ? numberLine("nmax_cal", "nmaxCal", "CAL n_max") : "",
    value("protocol") === "CAL" ? numberLine("u_cal", "uCal", "CAL u") : ""
  ]));
  const protocol = value("protocol");
  if (protocol) parts.push(section("protocol", [`name = ${tomlString(protocol)}`]));
  parts.push(section("requirements", [numberLine("target_key_bps", "targetKey", "target key")]));
  parts.push(`[reach]\nenabled = ${value("reach")}\n`);
  if (mode === "spectral") {
    const laserSpectrum = csvSpec("laser", "laser");
    parts.push(section("laser", [
      value("laserModel") ? `model = ${tomlString(value("laserModel"))}` : "",
      numberLine("lorentz_width_hz", "laserLinewidth", "laser linewidth"),
      laserSpectrum ? `spectrum = ${laserSpectrum}` : "",
      numberLine("r3", "r3", "r3"),
      numberLine("r2", "r2", "r2"),
      numberLine("fc_hz", "fc", "laser fc"),
      numberLine("C4", "c4", "C4"),
      numberLine("C3", "c3", "C3"),
      numberLine("C2", "c2", "C2"),
      numberLine("B_hz", "bHz", "B"),
      numberLine("gamma", "gamma", "gamma"),
      numberLine("delta", "delta", "delta")
    ]));
    const lineSpectrum = csvSpec("line", "line");
    if (lineSpectrum) lineLines.push(`spectrum = ${lineSpectrum}`);
    parts.push(section("scheme", [
      value("schemeLasers") ? `lasers = ${tomlString(value("schemeLasers"))}` : "",
      value("schemeCompensation") ? `compensation = ${tomlString(value("schemeCompensation"))}` : ""
    ]));
    parts.push(section("physics", [
      numberLine("lambda_s_nm", "lambdaS", "lambda_s"),
      numberLine("lambda_q_nm", "lambdaQ", "lambda_q"),
      numberLine("s0", "s0", "s0"),
      numberLine("fc2_hz", "fc2", "fc2"),
      numberLine("n", "indexN", "n"),
      numberLine("K", "kFactor", "K")
    ]));
    const actuatorLines = [];
    if (value("actuatorInput") === "omega") actuatorLines.push(numberLine("omega_a_rad_s", "omegaA", "omega_a"));
    if (value("actuatorInput") === "bandwidth") actuatorLines.push(numberLine("bandwidth_hz", "actuatorBandwidth", "actuator bandwidth"));
    if (value("actuatorInput") === "csv") {
      const response = actuatorResponseSpec();
      if (response) actuatorLines.push(`response = ${response}`);
    }
    parts.push(section("actuator", actuatorLines));
    parts.push(section("operation", [
      numberLine("sigma_limit_rad", "sigmaLimit", "sigma threshold"),
      numberLine("tau_max_s", "tauMax", "tau cap"),
      numberLine("tau_ps_s", "operationTauPs", "operation tau_PS")
    ]));
  }
  parts.splice(1, 0, section("line", lineLines));
  const extra = value("extraToml");
  if (extra) parts.push(extra);
  $("toml").value = parts.filter(Boolean).join("\n").trim() + "\n";
  builtText = $("toml").value; editorDirty = false;
}

async function readFiles(inputId) {
  const files = Array.from($(inputId).files || []);
  return Promise.all(files.map((file) => new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve({ name: file.name, content: String(reader.result) });
    reader.onerror = () => reject(reader.error);
    reader.readAsText(file);
  })));
}

function showStatus(text) {
  $("status").textContent = text;
}

function downloads(artifacts) {
  $("downloads").replaceChildren();
  for (const [name, href] of Object.entries(artifacts || {})) {
    const link = document.createElement("a");
    link.href = href;
    const names = {"report.html": "HTML report", "report.md": "Markdown report", "result.json": "Full results (JSON)", "spectrum.png": "Phase spectrum", "cumulative_variance.png": "Accumulated variance"};
    const parts = name.split("/");
    const prefix = parts.length > 1 ? "Variant " + parts[0].replace("variant_", "") + " · " : "";
    link.textContent = prefix + (names[parts.at(-1)] || parts.at(-1));
    link.download = name;
    $("downloads").appendChild(link);
  }
}

function renderReport(htmlHref) {
  if (htmlHref) {
    $("report").src = htmlHref;
    $("report").removeAttribute("srcdoc");
  } else {
    $("report").removeAttribute("src");
    $("report").srcdoc = "";
  }
}

async function postJSON(path, payload) {
  const response = await fetch(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload)
  });
  const data = await response.json();
  if (!response.ok || !data.ok) {
    showFieldErrors(data.field_errors);
    if (data.missing_inputs) showFieldErrors(Object.fromEntries(data.missing_inputs.map((path) => [path, "Required explicit input."])));
    if (editorDirty) $("configurationAdvanced").open = true;
    const lines = [data.error || "Calculation failed"];
    if (data.missing_inputs) lines.push("Missing inputs: " + data.missing_inputs.join(", "));
    if (data.diagnostics) lines.push(JSON.stringify(data.diagnostics, null, 2));
    throw new Error(lines.join("\n"));
  }
  return data;
}

async function runSingle() {
  if (!editorDirty) buildToml();
  showStatus("Running...");
  downloads({});
  renderReport(null);
  const files = await readFiles("files");
  const input = { toml: $("toml").value, files };
  if (!editorDirty && $("toml").value === builtText) {
    input.intensities_per_user = {};
    const controls = value("protocol") === "CAL" ? {u_cal:"uCal"} : {decoy_big:"decoyBig",decoy_medium:"decoyMedium",decoy_mini:"decoyMini"};
    for (const [key, id] of Object.entries(controls)) if (value(id)) input.intensities_per_user[key] = Number(value(id));
    if (exampleSource) input.example_source = exampleSource;
  }
  const data = await postJSON("/api/run", input);

  if (data.result.web_input) {
    for (const [path, notice] of Object.entries(data.result.web_input.warnings)) {
      for (const [id, spec] of Object.entries(fields)) if (spec.path === path) fieldMessage(id, notice, "warning");
    }
  }
  downloads(data.artifacts);
  renderReport(data.artifacts["report.html"]);
  const key = data.result && data.result.selected ? data.result.selected.key_bps : undefined;
  $("status").scrollIntoView({block: "start"});
  showStatus(`Complete: ${data.kind}\nKey rate: ${key === undefined ? "see report" : key} bit/s`);
}

async function runCompare() {
  showStatus("Running comparison...");
  downloads({});
  renderReport(null);
  const variants = [
    { label: value("labelA"), toml: $("tomlA").value, files: [] },
    { label: value("labelB"), toml: $("tomlB").value, files: [] }
  ];
  if (value("tomlC")) variants.push({ label: value("labelC"), toml: $("tomlC").value, files: [] });
  const items = value("lengths").split(",");
  const lengths = items.map((item) => Number(item.trim()));
  if (items.some((item) => !item.trim()) || lengths.some((n) => !Number.isFinite(n) || n <= 0)) {
    fieldMessage("lengths", "All working total lengths must be finite and > 0 [km].");
    throw new Error("Correct the working lengths before comparing.");
  }
  fieldMessage("lengths", "");
  const data = await postJSON("/api/compare", { variants, working_total_lengths_km: lengths, workers: 1 });
  downloads(data.artifacts);
  renderReport(data.artifacts["report.html"]);
  $("status").scrollIntoView({block: "start"});
  showStatus(`Complete: ${data.kind}\nVariants: ${variants.length}`);
}

document.querySelectorAll(".tab").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll(".tab").forEach((item) => item.classList.remove("active"));
    document.querySelectorAll(".panel").forEach((item) => item.classList.remove("active"));
    button.classList.add("active");
    $(button.dataset.panel).classList.add("active");
  });
});

["mode", "protocol", "detectorMode", "laserModel", "schemeCompensation"].forEach((id) => {
  $(id).addEventListener("change", updateMode);
});
$("build").addEventListener("click", () => {
  try {
    buildToml();
    showStatus("Configuration built.");
  } catch (error) {
    showStatus(error.message);
  }
});
$("run").addEventListener("click", () => runSingle().catch((error) => showStatus(error.message)));
$("compareRun").addEventListener("click", () => runCompare().catch((error) => showStatus(error.message)));

async function initialize() {
  $("run").disabled = true; $("build").disabled = true; $("loadExample").disabled = true;
  const response = await fetch("/api/schema"); if (!response.ok) throw new Error("Cannot load input validation schema.");
  const data = await response.json(); fields = data.fields; ruleMessages = data.rules;
  const tooltip = document.createElement("div"); tooltip.id = "fieldTooltip"; tooltip.className = "field-tooltip"; tooltip.hidden = true; tooltip.setAttribute("role", "tooltip"); document.body.append(tooltip);
  decorateFields();
  for (const id of Object.keys(fields)) {
    const el = $(id); if (!el) continue;
    el.addEventListener("input", () => {
      if (id === "toml") { editorDirty = true; return; }
      if (["labelA","labelB","labelC","tomlA","tomlB","tomlC","lengths","exampleChoice"].includes(id)) return;
      editorDirty = false;
      if (exampleSource) $("exampleSource").textContent = exampleSource + " Form values have been edited after loading; review the assumptions.";
      updateMode(); updateLossPriority(); validateForm(); updateLossPriority();
    });
  }
  $("loadExample").addEventListener("click", () => loadExample().catch((error) => showStatus(error.message)));
  $("run").disabled = false; $("build").disabled = false; $("loadExample").disabled = false;
  updateMode(); buildToml();
}
initialize().catch((error) => showStatus(error.message));
