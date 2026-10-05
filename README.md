<a id="english"></a>

# TF-QKD and mode-pairing noise model

[Русский](#russian) · [User guide](docs/USER_GUIDE.md) · [Validation](docs/VALIDATION.md) · [Sources](SOURCES.md)

Estimate laser and fiber phase noise in quantum key distribution, compare
stabilization schemes, and calculate asymptotic key rates from specified or
measured inputs.

**TF-QKD:** calculate phase-noise spectra, phase variance, the usable transmission
window, duty cycle, and SNS-AOPP or CAL key rate. Compare uncompensated fiber,
classical round-trip compensation, dual-band stabilization, and ideal fiber-noise
cancellation. Measured spectra and residual phase RMS are also supported.

**Mode-pairing QKD:** calculate phase increments between paired pulses and a
base asymptotic key bound from explicit detection gains. The phase calculation
and key-rate calculation are separate commands. A key estimate based on observed
gains is conditional on those observations.

## Get started

Requires Python **3.11+**. Run these commands in a terminal:

```bash
git clone https://github.com/Quarzii/tfqkd-lab.git
cd tfqkd-lab
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m tfqkd.lab_web --port 8765
```

On Windows, activate the environment with `.venv\Scripts\activate`.
Open <http://127.0.0.1:8765> in your browser. Stop the server with Ctrl+C.
The interface runs locally and offers TF-QKD input forms and reports.

For a TF-QKD calculation from a file:

```bash
python -m tfqkd.lab_run examples/laboratory/stabilized_dual.toml --output results/my_run
```

This writes Markdown, HTML, JSON, and, for spectral calculations, plots.
Replace example parameters with your installation's values; example coefficients
are calculation inputs, not equipment defaults.

For mode-pairing phase noise and key rate:

```bash
python -m tfqkd.mp_phase examples/mode_pairing/phase.toml --output results/my_mp_phase.json
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
```

## Choose your inputs

- **Noise prediction:** supply explicit laser/fiber coefficients or measured PSD files.
- **Measured phase:** supply residual phase RMS and its operating window to estimate the key rate at that measured condition.
- **Mode-pairing key rate:** supply detection gains, error gains, pairing parameters, and phase or drift inputs.

[The user guide](docs/USER_GUIDE.md) explains units, outputs, and how to select an
example. [Input reference](docs/UNIVERSAL_INPUTS.md) covers TF-QKD configuration
fields and CSV formats.

## Scope and accuracy

The TF-QKD model follows Bertaina et al.; classical compensation uses the
Williams round-trip model. The mode-pairing calculation uses the base Zeng
protocol and documented phase/decoy assumptions. See [sources](SOURCES.md).

Key rates are asymptotic: finite-size security analysis is outside the model.
Classical compensation omits phase-measurement detection noise, so its rate is
optimistic. Results depend on the supplied spectra, timing, and apparatus values;
coefficients from one route are not calibrated predictions for another.
Missing noise amplitudes or intrinsic-error parameters can produce a labelled
upper model estimate rather than an installation forecast.

Numerical consistency checks and published-experiment comparisons have different
scopes. The tool does not reproduce every published rate or measured spectrum.
See [validation and discrepancies](docs/VALIDATION.md) and
[model limitations](docs/UNKNOWNS.md).

## Check the installation

```bash
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
python scripts/validate/validate_stage3_a.py --calibrate
```

These checks run without article PDFs. Paper-dependent reproduction checks are
explained in [SOURCES.md](SOURCES.md).

## License

[MIT](LICENSE) covers original contributions. Reference software, adapted author
functions, and measurement files retain their own licenses and attribution:
[third-party notices](THIRD_PARTY_NOTICES.md).

---

<a id="russian"></a>

# Модель шумов TF-QKD и mode-pairing QKD

[English](#english) · [Руководство](docs/USER_GUIDE.md) · [Проверки модели](docs/VALIDATION.md) · [Источники](SOURCES.md)

Модель оценивает фазовый шум лазеров и оптоволокна в квантовом распределении
ключей, сравнивает способы стабилизации и рассчитывает асимптотическую скорость
ключа по заданным параметрам или измерениям.

**TF-QKD:** спектры и дисперсия фазового шума, рабочее окно передачи, доля
времени передачи ключа, скорость SNS-AOPP или CAL. Можно сравнить волокно без
компенсации, классическую компенсацию с двойным проходом, двухполосную
стабилизацию и идеальное устранение шума волокна. Поддерживаются измеренные
спектры и остаточная СКО фазы.

**Mode-pairing QKD:** приращения фазы между спаренными импульсами и базовая
асимптотическая оценка скорости ключа по явно заданным вероятностям регистрации.
Расчёт фазы и расчёт ключа запускаются отдельно. Оценка по наблюдаемым
регистрациям обусловлена этими наблюдениями.

## Запуск

Нужен Python **3.11+**. Выполните в терминале:

```bash
git clone https://github.com/Quarzii/tfqkd-lab.git
cd tfqkd-lab
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m tfqkd.lab_web --port 8765
```

В Windows активируйте окружение командой `.venv\Scripts\activate`.
Откройте <http://127.0.0.1:8765> в браузере. Сервер останавливается сочетанием
Ctrl+C. Локальный интерфейс позволяет задать входные данные TF-QKD и получить отчёт.

Расчёт TF-QKD из файла:

```bash
python -m tfqkd.lab_run examples/laboratory/stabilized_dual.toml --output results/my_run
```

Результат сохраняется в Markdown, HTML и JSON; для спектрального расчёта также
строятся графики. Замените параметры примера своими: опубликованные коэффициенты
не являются значениями по умолчанию для вашей установки.

Расчёт фазового шума и скорости ключа mode-pairing:

```bash
python -m tfqkd.mp_phase examples/mode_pairing/phase.toml --output results/my_mp_phase.json
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
```

## Какие данные нужны

- **Прогноз шума:** коэффициенты модели лазеров и волокна либо измеренные спектры PSD.
- **Оценка по измеренной фазе:** остаточная СКО фазы и соответствующее рабочее окно.
- **Скорость mode-pairing:** вероятности регистраций и ошибок, параметры спаривания, данные о фазе или дрейфе.

[Руководство](docs/USER_GUIDE.md) объясняет выбор примера, единицы и результаты.
[Справочник входных данных](docs/UNIVERSAL_INPUTS.md) описывает поля TF-QKD и
форматы CSV. Подробная документация — на английском.

## Область применения

Расчёт TF-QKD основан на модели Bertaina и соавторов, классическая компенсация —
на модели Williams. Для mode-pairing используется базовый протокол Zeng с
описанными допущениями о фазе и decoy-оценках. См. [источники](SOURCES.md).

Скорости ключа асимптотические: анализ безопасности конечной сессии не включён.
В классической компенсации отсутствует шум измерения фазы, поэтому её оценка
оптимистична. Прогноз зависит от спектров, временных параметров и данных установки;
коэффициенты одной трассы не дают калиброванный прогноз для другой.
При неизвестных амплитудах шума или собственных ошибках результат может быть
помечен как верхняя оценка модели.

Численная согласованность и сравнение с экспериментом проверяют разные свойства.
Модель воспроизводит не все опубликованные скорости и формы спектров.
[Проверки и расхождения](docs/VALIDATION.md) и
[ограничения](docs/UNKNOWNS.md) описаны отдельно.

## Проверка установки

```bash
python -m unittest discover -s tests
TFQKD_FULL_TESTS=1 python -m unittest discover -s tests
python scripts/validate/validate_stage3_a.py --calibrate
```

Эти проверки не требуют PDF статей. Для проверок по исходным публикациям
см. [SOURCES.md](SOURCES.md).

## Лицензия

[MIT](LICENSE) распространяется на оригинальные части проекта. Авторский код,
его адаптации и измерительные данные сохраняют свои лицензии и атрибуцию:
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
