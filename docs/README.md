# Twin-Field QKD: опубликованная модель и компенсация Williams

## Ввод измеренной фазы и спектров

[Формат новых входов](MEASURED_INPUTS.md): остаточная СКО фазы с окном
измерения, CSV без подгонки, полные потери плеч и два детектора отдельно.
Проверенный отчёт (`../reports/archive/MEASURED_INPUTS_REPORT.md`; research archive not included in the public snapshot) содержит условное сравнение
Pittaluga/Zhou и численные расхождения; прежняя привязка СКО Zhou к трём
длинам отозвана. Физическое ядро сохранено. Часть D отложена.

```bash
python -m tfqkd.lab_run examples/measured_inputs/pittaluga2021_368.702_phase_bound.toml --output results/my_phase_run
python scripts/validate/validate_measured_inputs.py --quick --regression --experiments
python scripts/validate/validate_measured_inputs.py --full
python scripts/report/report_measured_inputs.py
```

## Аудит

Отчёт аудита (`../reports/archive/AUDIT_REPORT.md`; research archive not included in the public snapshot): независимая реализация без импорта ядра,
извлечение опубликованных входов Pittaluga2021/Clivati2022/Zhou2023 и ревизия
классического блока. Сквозная экспериментальная валидация скорости пока не
установлена; неполные входы и причины отказа сохранены в
examples/audit (`../examples/audit/README.md`; research archive not included in the public snapshot). Часть D отложена.

```bash
python scripts/validate/audit_independent.py
python scripts/validate/audit_experiments.py
python scripts/run/run_audit.py
python scripts/report/report_audit.py
```

## Универсальный ввод установки

Библиотека пресетов исключена из инструмента. Параметры задаются явно,
паспортными величинами или измерительным CSV (PSD лазера/трассы;
частотная характеристика исполнителя). Отсутствующие амплитуды шума дают
верхнюю оценку скорости и отдельные условные требования к оборудованию.
Другие обязательные параметры не заполняются автоматически.

```bash
python -m tfqkd.lab_run examples/bertaina2024_table3.toml --output results/my_run
python -m tfqkd.lab_inputs my_setup.toml
python -m tfqkd.measured_inputs examples/b6/jiang2008_urban86.json --node line
python scripts/validate/validate_universal_inputs.py --quick --inputs --regression
```

[Формат ввода и примеры](UNIVERSAL_INPUTS.md).
Проверенный отчёт доработки B (`../reports/archive/UNIVERSAL_INPUTS_REPORT.md`; research archive not included in the public snapshot).
`examples/bertaina2024_table3.toml` — явно заданная опубликованная установка,
**не значения по умолчанию**. `configs/config.toml` используется только старыми
расчётами воспроизведения источников и их тестами; новый resolver не берёт
из него физические параметры. `n` имеет отдельно документированное справочное
умолчание из паспорта Corning. Часть C добавляет ранжирование чувствительности, сравнение вариантов,
разложение дисперсии и отчёты Markdown/HTML; интерфейс относится к части D.

Исторические отчёт B (`../reports/archive/STAGE3_B_REPORT.md`; research archive not included in the public snapshot),
доработки B6 (`../reports/archive/STAGE3_B_AMENDMENTS_REPORT.md`; research archive not included in the public snapshot) и аудит источников (`../reports/archive/B6_SOURCE_AUDIT.md`; research archive not included in the public snapshot)
сохранены; их прежние команды с выбором пресетов заменены явным вводом.

Python 3.11+, numpy, scipy, matplotlib. Статьи в sources/papers/, данные в sources/data/.
Все параметры исторических проверок и единицы заданы в configs/config.toml; расчёты сохраняют снимок входов.
Используются предоставленные открытые источники. Внутренние параметры установки
и новые лабораторные измерения от пользователя не требуются.

## Расчёт сценариев и графиков на выбранной сетке

```bash
python scripts/run/run_stage1.py --config configs/config.toml
python -m tfqkd.report
python scripts/run/run_keyrates.py --config configs/config.toml
python scripts/validate/validate_keyrates.py
python scripts/run/run_stage2_checks.py --config configs/config.toml
python scripts/run/run_classical.py --config configs/config.toml
python -m unittest discover -s tests -v > results/all_tests.txt 2>&1
python scripts/report/report_classical.py
```

`scripts/run/run_classical.py` запускает актуальный `scripts/run/run_classical_actuator.py`: исполнительный
полюс, решаемое время передачи и ноль скорости ключа. Длительный скан сохраняет
промежуточный actuator_checkpoint.json и сообщает прогресс по завершении кривых.
Отчёты: STAGE1_REPORT.md, STAGE2_REPORT.md, BLOCK5_REPORT.md.
Т4 содержит количественное расхождение с измерениями и остаётся пределом валидации.

## Текущий блок 5

- C(s)=g/[s*(1+s/omega_a)] — явное инженерное допущение пользователя.
- G=C(s)*(1+exp(-s*tau_RT)) — Williams Eq. A6; удалённый PSD — Eqs. A3/A8.
- g_crit(omega_a,L) определяется численно по мнимому корню и проверяется по исходной комплексной петле.
- Неустойчивый/предельный g выдаёт ошибку. В пакетной таблице скорость пустая с явным статусом.
- tau_Q находится из Eq. 4; 100 ms — верхняя отсечка, а не фиксированное рабочее время.
- Рабочий цикл, SNS-AOPP/CAL и потери рассчитываются из найденного времени и фактической sigma.
- Нули подписанной нижней оценки ключа определяют дальность. Устойчивость учитывается отдельно.
- При фиксированном g/g_crit абсолютное g меняется с L; дополнительно исследовано фиксированное g.

Графики: results/actuator_stability.png, actuator_suppression.png,
actuator_classical_vs_dual.png, actuator_key_reach.png.
Данные: actuator_gain_scan.csv, actuator_keyrate_curves.csv, actuator_dual_keyrate.csv,
actuator_key_reach.csv, actuator_operating_convergence.csv, actuator_spectra.npz.
Полные параметры и сводные проверки: results/actuator_evidence.json.

## Отдельные модули

```bash
python -m tfqkd.config
python -m tfqkd.transfers
python -m tfqkd.spectra
python -m tfqkd.integration
python -m tfqkd.protocol
python -m tfqkd.noise
python -m tfqkd.reference
python -m tfqkd.keyrates
python -m tfqkd.linewidth
python -m tfqkd.classical
python -m tfqkd.classical_keyrate
```

Модули принимают `--config`. classical показывает также формальную идеальную
петлю; classical_keyrate выводит актуальные рабочие точки с исполнительным полюсом.

## Исторические расчёты

results/stage2_fixed_time/ содержит полный прежний расчёт с C=g/s и фиксированными
100 ms. Он отвечает на вопрос о фазовом пределе за заданное время и не является
актуальной дальностью QKD. Его диагностический скрипт — scripts/run/run_classical_fixed_time.py.
Отсутствие конечного g_crit у идеального интегратора — вырождение этой модели;
его аналитическое объяснение сохранено в текущем отчёте.
Другие архивы: results/stage1_original/ и results/stage2_before_integrator/.
Ограничения источников: UNKNOWNS.md и PUBLIC_SOURCE_LIMITS.md.

## Этап 3, часть A: ускоренное ядро

По умолчанию `grid.mode="fast"`: базовая сетка **4097** точек. Режим
`reference` сохраняет **65537** точек. В классической схеме к базовой сетке
добавляются точки около найденных полюсов; окончательное число точек зависит
от параметров и проверяется удвоением сетки. Неустойчивость или отсутствие
сходимости приводят к явной ошибке.

```bash
python -m tfqkd.integration --grid-mode fast
python -m tfqkd.classical_keyrate --grid-mode fast
python -m tfqkd.engine --grid-mode fast
python scripts/validate/validate_stage3_a.py --quick
python scripts/validate/validate_stage3_a.py --calibrate
python scripts/validate/validate_stage3_a.py --full
python scripts/run/benchmark_stage3_a_before.py
python scripts/run/benchmark_stage3_a.py
python scripts/report/report_stage3_a.py
```

`--quick` выполняет быстрые юнит-тесты. `--full` воспроизводит Т1–Т7,
Figure 3 и авторские протоколы на reference-сетке, сохраняя результаты в
`results/stage3_a/full/`. Принятое ограничение Т4 сохраняется.
`--calibrate` сверяет вложенные сетки для всех семи сценариев Table I.
Пакетное сканирование: `tfqkd.engine.scan(config)`; выходные массивы имеют
оси `(omega_a, g, L)`. Спектральные рабочие массивы разбиты на ограниченные
по памяти блоки произведения g×L с последней осью частоты. Процессы используются
для независимых вариантов omega_a; внутри блока работают операции NumPy.

`performance.workers` задаёт число процессов. Для запуска из Python требуется
обычный `if __name__ == '__main__':` при создании пула процессов; это условие
выполнено в предоставленных скриптах. `workers=1` запускает тот же расчёт
последовательно. `batch_memory_mb` задаёт оценочный бюджет рабочих массивов
блока первой проверки; это не жёсткий лимит RSS при дальнейших уточнениях.
Кэш ограничен числом записей `cache_entries`; его можно очистить через
`tfqkd.cache.clear_spectral_cache()`.

4097 — минимальная прошедшая сетка **среди проверенных вложенных 2^k+1**
для зафиксированных сценариев и времён. Это не доказательство для произвольных
пользовательских спектров. Для проверки новых режимов доступны reference и
повторная калибровка.

## Этап 3, часть B: обновлённый ввод лаборатории

[UNIVERSAL_INPUTS.md](UNIVERSAL_INPUTS.md) описывает прямые значения, паспортные
пересчёты, CSV, условные допуски и четырёхсторонний потолок выигрыша.
`configs/lab_defaults.toml` содержит вычислительные настройки и разрешённые правила
поиска, а не параметры установки. Основная команда — `python -m tfqkd.lab_run`.
Исторические проверки: `python scripts/validate/validate_universal_inputs.py --full`.
Документация части C: [PART_C_OUTPUTS.md](PART_C_OUTPUTS.md).
Проверенный отчёт части C (`../reports/archive/STAGE3_C_REPORT.md`; research archive not included in the public snapshot) содержит результаты тестов,
замеры и ссылки на Markdown/HTML примеры одного запуска и сравнения.
Часть D начинается после приёмки части C.
