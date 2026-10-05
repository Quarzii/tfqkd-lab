"""Render the report solely from saved execution evidence: python -m tfqkd.report."""

import json
from .config import ROOT


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    e = json.loads((ROOT / "results/evidence.json").read_text())
    t1, t2, c = e["T1"], e["T2"], e["config"]
    scenarios = table(["Сценарий", "tau_Q, с", "Авторы, с", "Разность, %", "sigma_phi, рад", "e_phi", "d", "Ограничение"],
                      [[r["scenario"], f'{r["tau_s"]:.9g}', f'{r["author_tau_capped_s"]:.9g}',
                        f'{r["tau_difference_percent"]:.6g}', f'{r["sigma_rad"]:.9g}', f'{r["e_phi"]:.9g}',
                        f'{r["duty"]:.9g}', r["status"]] for r in e["tableI"]])
    fiber = table(["Диапазон, Гц [a,b)", "N", "Медиана измерение/модель", "Уровень, дБ", "RMS log10", "Наклон измерения", "Наклон модели", "Разность наклонов"],
                  [[f'{r["low_hz"]:g}–{r["high_hz_exclusive"]:g}', r["points"], f'{r["median_measurement_over_model"]:.7g}',
                    f'{r["median_residual_db"]:.6g}', f'{r["rms_log10_residual"]:.6g}', f'{r["measured_slope"]:.6g}',
                    f'{r["model_slope"]:.6g}', f'{r["slope_difference"]:.6g}'] for r in e["T4"]])
    refinement = table(["Точек авторской сетки", "Авторское tau, с", "Разность с уточнённым интегралом, %"],
                       [[r["points"], f'{r["author_tau_s"]:.10g}', f'{r["difference_percent"]:.7g}']
                        for r in e["figure2"]["author_refinement_at_worst_point"]])
    finite = t2["finite_record_expectations_eq38"]
    variance_table = table(["Оценка", "Среднее", "Эмпирическая SE", "Разность с интегралом, %", "Разность / SE"],
                          [[key, f'{t2[key]["mean"]:.12g}', f'{t2[key]["standard_error"]:.8g}',
                            f'{t2[key]["difference_percent"]:.7g}', f'{t2[key]["difference_in_standard_errors"]:.6g}']
                           for key in ("zero_mean_second_moment", "centered_variance")])
    max_grid = max(r["variance_max_relative_change"] for r in e["convergence"])
    max_spectrum = max(r["max_psd_relative_error"] for r in e["tableI"])
    max_fig8 = max(abs(r["difference_percent"]) for r in e["figure8"])
    text = f'''# Отчёт этапа 1

Реализованы блоки 1–3. Выполнены Т1, Т2, части Т3 для Table I / Figures 2 и 8, а также Т4.
**Этап 1 принят пользователем с уточнёнными критериями.** Т2 теперь сравнивается с интегралом от 1/T и конечным ожиданием. Строгая область Т4 ограничена f<100 Гц; численные невязки сохранены.
Figure 3 перенесена пользователем на этап 2. Блоки 4–5 и Т5–Т7 не выполнялись.
Скрипт возвращает 0 при успешных Т1, уточнённом Т2 и тестах кода; это не утверждение точного экспериментального согласия Т4.

## Воспроизводимость

Команды: `python run_stage1.py --config config.toml`, `python -m tfqkd.report`.
Первый скрипт также выполняет unittest. Версии: {e['versions']}.
Снимок входов, статусы и SHA-256 источников: [evidence.json](results/evidence.json).
Эталон — исходные функции ячейки 3 QKD.ipynb, извлечённые через AST без изменения тел.
Графики построены по вычисленным массивам; пиксели печатных рисунков не оцифровывались.

## Реализованные уравнения

| Модуль | Формула | Источник |
|---|---|---|
| spectra | Свободный лазер | Eq. F1, bertaina2024 |
| spectra | Стабилизированный лазер | Eq. F2, bertaina2024 |
| spectra | Шум резонатора | Eq. F3, bertaina2024 |
| transfers | G(f), G0 | Eq. F4 и следующий абзац, bertaina2024; QKD.ipynb G1 |
| spectra | Свободное волокно | Eq. 6, bertaina2024 |
| spectra | Подавленный шум волокна | Eq. 8, bertaina2024; QKD.ipynb stabfibnoise |
| spectra | Шум детектирования | Appendix G, ненумерованная S_detection, bertaina2024; QKD.ipynb fiberdetection |
| transfers, spectra | Общий / независимые лазеры | Eqs. 5, 7, bertaina2024; QKD.ipynb calc_spectra |
| integration | Дисперсия и накопленный интеграл | Eq. 4, bertaina2024; конечный верхний предел по заданию |
| integration | Порог времени и ограничение | Sec. V, bertaina2024; QKD.ipynb calc_sigma_tau / plot_panel_scenarios |
| protocol | e_phi = sigma_phi²/4 | Eq. 1, bertaina2024, запрошенное приближение |
| protocol | Рабочий цикл | Sec. I, p. 3, ненумерованная формула; QKD.ipynb duty_stabilization |
| noise | FIR и линейная свёртка | Eqs. 104, 37, Sec. III-B, VI-C, Appendix II, kasdin1995 |
| noise | Дискретный PSD и степенной предел | Eqs. 98, 99, kasdin1995; двусторонняя конвенция Eq. 15 |
| noise | Стационарная и конечная дисперсии | Eq. 111 и вывод из Eqs. 37–38, kasdin1995 |

Ненумерованным формулам не присвоены вымышленные номера. Коэффициенты — Table III / QKD.ipynb.
G0 вычисляется как в G1, без замены округлённым табличным 3.55·10¹³.
n={c['physics']['n']}, K={c['physics']['K']}. В Eq. 7 коэффициент волокна равен 1;
K применяется к Eq. 5. Шум детектирования добавлен один раз без K при стабилизации волокна.

## Т1: нормировка

alpha={c['validation']['alpha']}, Qd={c['validation']['input_variance']}, fs={c['validation']['sample_rate_hz']:g} Гц,
N={c['validation']['samples']}, реализаций={c['validation']['realizations']}, seed={c['validation']['seed']}.
Welch: {c['validation']['welch_window']}, сегмент {c['validation']['welch_segment']}, перекрытие {c['validation']['welch_overlap']},
`detrend=False`, `scaling="density"`, `return_onesided=True`.
Использована линейная свёртка с нулевым прошлым без последующей перенормировки ряда.

Положительная часть двусторонней Eq. 98 удвоена; бин Найквиста не удвоен, DC исключён.
В полосе {c['validation']['comparison_min_hz']:g}–{c['validation']['comparison_max_hz']:g} Гц:

- Welch / Eq. 98 = **{t1['discrete_psd_ratio']['mean']:.9f} ± {t1['discrete_psd_ratio']['standard_error']:.9f} SE**.
- Welch / Eq. 99 = **{t1['power_law_ratio']['mean']:.9f} ± {t1['power_law_ratio']['standard_error']:.9f} SE**.
- Критерий ±{c['validation']['confidence_sigma']:g} эмпирических SE для Eq. 98: **{t1['pass']}**.

SE получена по независимым реализациям. Eq. 98 — точный дискретный спектр алгоритма;
Eq. 99 — его низкочастотное приближение (kasdin1995, Sec. VI-B, p. 820).
[График Т1](results/T1_psd.png), [данные](results/T1_psd.csv).

## Т2: дисперсия

Интеграл односторонней Eq. 98 от 0 до fs/2: **{t2['psd_integral']:.12g}**,
оценка ошибки квадратуры {t2['quadrature_error']:.3g}.
Стационарная Eq. 111: **{t2['stationary_variance_eq111']:.12g}**.
Единицы синтетические, не экспериментальные rad².

{variance_table}

**Прежний Т2 (диагностика):** сравнение со стационарным полным интегралом не прошло.

Уточнённый Т2: T={t2["revised"]["record_time_s"]} с, 1/T={t2["revised"]["lower_hz"]} Гц.
Интеграл от 1/T до Найквиста = **{t2["revised"]["band_integral"]:.12g}**,
разность с измеренной центрированной дисперсией = **{t2["revised"]["difference_in_standard_errors"]:.6g} SE**.
Вместе с конечным ожиданием оба критерия 3 SE: **{t2["revised"]["pass"]}**.
Жёсткий срез 1/T — заданный критерий; вычитание среднего не объявляется идеальным частотным фильтром.
Оценка mean(x²) при заданном нулевом ансамблевом среднем (Kasdin Definition 1 / Eq. 41)
находится в пределах 3 SE, но это не заменяет проваленную проверку центрированной оценки.

Диагностика конечного ряда по Eqs. 37–38, kasdin1995:
E[mean(x²)]={finite['mean_square']:.12g}, Var(mean(x))={finite['variance_of_record_mean']:.12g},
E[np.var(x)]={finite['centered_variance']:.12g}.
Разность наблюдаемой центрированной дисперсии и её конечного ожидания:
**{t2['finite_centered_difference_in_standard_errors']:.6g} SE**.
Эти ожидания получены из весов того же FIR-фильтра; ряд и параметры не менялись после проверки.
[Реализации Т2](results/T2_realizations.csv).

## Т3: Table I и Figures 2, 8

Table I содержит конфигурации, не табулированные времена. Колонка «Авторы» — запуск calc_sigma_tau
на его сетке {c['grid']['author_points']} точек с ограничением времени.
LB={c['operation']['LB_km']:g} км, sigma_limit={c['operation']['sigma_limit_rad']:g} рад,
tau_max={c['operation']['tau_max_s']:g} с, tau_PS={c['operation']['tau_ps_s']:g} с.

{scenarios}

Максимум относительного расхождения PSD с авторами по семи сценариям: **{max_spectrum:.6g}**.
Для Figure 8 сценарий 3 взят со стабилизацией волокна, как в ноутбуке;
при ANY в сценариях 6–7 использовано авторское dL=2.5 км.

Figure 8: максимум модуля разности sigma² с авторской сеткой 1000 точек на {len(e['figure8'])} контрольных парах
(времена {c['validation']['comparison_times_s']} с) — **{max_fig8:.6g}%**.
Правая ось авторского графика и код показывают sigma в rad; на графике воспроизведена sigma,
в численном сравнении используется sigma².
[Figure 8](results/figure8.png), [сравнение](results/figure8_comparison.csv).

Figure 2: {e['figure2']['points']} точек контуров, восемь комбинаций.
Максимальное расхождение времени с авторской сеткой 1000 точек — **{e['figure2']['max_abs_tau_difference_percent']:.7g}%**
при dL={e['figure2']['worst_author_difference']['delta_L_km']:.9g} км.
Уточнённый интеграл даёт {e['figure2']['worst_author_difference']['refined_tau_s']:.10g} с.
Проверка этой точки исходным авторским интегратором на разных сетках:

{refinement}

В calc_sigma применены правые прямоугольники, причём включён интервал непосредственно ниже
отображаемого нижнего предела. Здесь применены трапеции с частичным интервалом точно от 1/tau.
Изменялись только сетки, физические параметры сохранены. Это результат о дискретизации, не провал модели.
[Интеграл при фиксированном времени против плотности сетки](results/figure2_integral_convergence.png), [данные](results/figure2_integral_convergence.csv).
[Figure 2](results/figure2.png), [сравнение](results/figure2_comparison.csv).
Figure 3 перенесена пользователем на этап 2 и не выполнялась.

## Сходимость и верхняя граница

Логарифмическая сетка: {c['grid']['f_min_hz']:g}–{c['grid']['f_max_hz']:g} Гц,
{c['grid']['points']} точек, уточнение до {c['grid']['refined_points']}.
Максимальное относительное изменение sigma² на контрольных временах: **{max_grid:.9g}**.
Максимальное изменение tau на контурах Figure 2: **{e['figure2']['max_abs_grid_change_percent']:.9g}%**.
Оба значения меньше заданного относительного допуска {c['grid']['relative_tolerance']:g}.

1 МГц воспроизводит конечный предел calc_sigma, не бесконечный интеграл.
Дополнительные границы {c['validation']['extended_cutoffs_hz']} Гц при неизменном исходном tau
увеличили sigma² максимум на **{max(r['change_from_1MHz_percent'] for r in e['upper_cutoff_sensitivity']):.7g}%**.
[Зависимость от верхнего предела](results/upper_cutoff_sensitivity.csv).
[Накопленная дисперсия](results/cumulative_variance.png), [численные данные](results/cumulative_variance.csv).

## Т4: линия 114 км

l=44 rad² Hz km⁻¹, fc1=100 Гц, L=114 км, множитель двойного прохода 4.
Источники: bertaina2024 Eq. 6, Appendix G, Figure 6; QKD.ipynb plot_laser_fiber_spectra.
Именно функция авторов и Appendix G связывают TXT с измерением 114 км;
сам заголовок TXT длину не содержит. Установка и измерения: clivati2022, основной PDF,
Figure 2 / Figure 3c и p. 4; второй PDF — исправление сведений о финансировании.

Сравнение выполнено в исходных частотных точках. Наклоны — описательная регрессия log10(PSD)
на log10(f), без изменения коэффициентов модели. Диапазоны полуоткрытые, точка 1 МГц исключена из статистик.

{fiber}

**Т4: строгая область теперь только f<fc1=100 Гц.** Последняя строка охватывает все доступные точки ниже fc1. Невязка полной Eq. 6 сохраняется; частотный множитель нельзя отбросить во всём диапазоне f<fc1. Асимптота 1/f² справедлива при f≪fc1.
Выше fc1 остатки — диагностика. Возможные, но не установленные для этих данных причины: механические/акустические резонансы (bertaina2024 Sec. III-A) и фон измерения (Appendix G описывает detection floor стабилизированного спектра). Данные не позволяют приписать весь хвост конкретному механизму. Статистическая значимость остатков неизвестна: нет неопределённостей измерения.
[График Т4](results/T4_fiber114km.png), [остатки](results/T4_residuals.csv).

## Проверки кода

```text
{(ROOT / 'results/unit_tests.txt').read_text().strip()}
```

## Решения, не продиктованные источниками

1. Чтение источников из корня вместо отсутствующей sources/; материалы не перемещались.
2. Трапеции, частичный граничный интервал и поиск корня по log(tau) — численная реализация Eq. 4.
   Отличие от авторских прямоугольников раскрыто и проверено сгущением сеток.
3. Сетки, допуск, seed, число реализаций, Welch и полоса сравнения — численные настройки, не физические коэффициенты;
   все собраны в config.toml. Синтетический alpha=0.5 выбран в стационарной области alpha<1 (Kasdin Eq. 111),
   Qd=1 задаёт условную единицу теста, а не измеренный уровень шума установки.
4. ±3 эмпирических SE — инженерный критерий Monte Carlo, не формула Kasdin и не доверительная полоса каждого бина.
5. Конечные ожидания выведены из Eqs. 37–38; служат вторым критерием уточнённого Т2; прежнее сравнение сохранено.
6. Медианы, RMS и регрессии наклонов описывают невязки Т4; показаны все декады авторского диапазона.
   Новая физическая модель высокочастотного фона не вводилась.
7. Расширение верхнего предела — отдельная численная проверка, не изменение эталона 1 МГц.
8. e_phi рассчитано в приближении из задания. Точная гауссовская функция из авторских протоколов здесь не подменяет Eq. 1.

## Содержимое UNKNOWNS.md

{(ROOT / 'docs/UNKNOWNS.md').read_text()}
'''
    (ROOT / "reports/archive/STAGE1_REPORT.md").write_text(text)
    print("Wrote", ROOT / "reports/archive/STAGE1_REPORT.md")


if __name__ == "__main__":
    main()
