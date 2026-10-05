# Public snapshot scope

This repository is a reviewed export of the v1.1 research workspace. The model,
unit tests and numerical configuration are byte-for-byte unchanged. Documentation
links point to public papers or explicitly mark omitted research material.
Local absolute paths in retained evidence are replaced by relative paths, without
changing numerical results. The public repository starts with a fresh history;
old local author identities and deleted research files are not carried over.

## Included

- Python model and local HTML/CSS/JavaScript interface; unit tests and numerical configuration.
- Complete published/example TOML inputs and small CSV/JSON fitting fixtures.
- One synthetic reference PSD for the direct-input example.
- User documentation, bilingual project descriptions and supervisor summary.
- Selected small validation Markdown/JSON evidence, including the discrepancies.
- Bertaina's notebook and four measurement files under attributed CC BY 4.0.
- A small set of validation/reproduction scripts; paper-dependent scripts require separately downloaded sources.

## Excluded

- Article PDFs, copies of article figures and extracted full article texts.
- Large generated arrays, all scan-by-scan reports and most historical stage artifacts.
- Unrelated third-party MATLAB code and obsolete equipment presets.
- Run logs, environment directories, caches, editor locks and original Git history.
- Credentials and local installation data; supplied inputs come from public sources or labelled numerical examples.

The MIT licence is limited to original contributions. Author-code adaptations
and reference data retain their own attribution/licence; see
[third-party notices](../THIRD_PARTY_NOTICES.md).

Quick/full unit suites, a complete example run, and fast/reference calibration
were executed in the export. Historical source-based T1–T7 checks are documented
but have not been re-run here without the excluded PDFs.
