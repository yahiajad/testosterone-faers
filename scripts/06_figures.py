"""Figures 1–4 (PNG 300 dpi + PDF) into $FAERS_DIR/figures."""
import os, glob, duckdb, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = os.environ.get("FAERS_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))  # data root (set FAERS_DIR to override)
HERE = os.path.dirname(os.path.abspath(__file__)); F = f"{B}/figures"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
INK, SA, PC = "#222222", "#2b6cb0", "#c05621"
save = lambda fig, n: [fig.savefig(f"{F}/{n}.{x}", dpi=300, bbox_inches="tight") for x in ("png", "pdf")]

# Fig 1: data flow (counts computed from the database built by 03_analysis.py)
import importlib.util
spec = importlib.util.spec_from_file_location("a", f"{HERE}/03_analysis.py"); a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
a.flag("t_ps", a.TESTO, ["PS"], a.TESTO_EXCL)
q = lambda sql: a.c.execute(sql).fetchone()[0]
nq = len(glob.glob(f"{a.P}/demo/*.parquet"))
cases = {ev: q(f"""select count(*) from base where sex='M' and primaryid in (select * from t_ps) and primaryid in
              (select primaryid from reac where pt in ({",".join(f"'{p}'" for p in pts)}))""") for ev, pts in a.EVENTS.items()}
n_raw = q(f"select count(*) from read_parquet('{a.P}/demo/*.parquet')")
n_demo = q("select count(*) from demo"); n_base = q("select count(*) from base")
n_men = q("select count(*) from base where sex='M'")
n_t = q("select count(*) from base where sex='M' and primaryid in (select * from t_ps)")
steps = [(f"FAERS quarterly files ({nq} quarters)", f"{n_raw:,} report records"),
         ("Latest version per case; FDA-deleted cases removed", f"{n_demo:,} cases"),
         ("Exact duplicates collapsed\n(event date, age, sex, country, drugs, PTs)", f"{n_base:,} unique reports"),
         ("Male patients", f"{n_men:,} reports"),
         ("Testosterone as primary suspect", f"{n_t:,} reports"),
         ("Cases: sleep apnoea | polycythaemia", f"{cases['sleep_apnoea']:,} | {cases['polycythaemia']:,}")]
excl = [f"Superseded case versions and\nFDA-deleted cases: {n_raw - n_demo:,}", f"Exact duplicates: {n_demo - n_base:,}",
        f"Female or sex not\nrecorded: {n_base - n_men:,}", f"Testosterone not primary\nsuspect: {n_men - n_t:,}", None]
fig, ax = plt.subplots(figsize=(7.4, 6.2)); ax.axis("off")
for i, (t, n) in enumerate(steps):
    y = 1 - i * 0.18
    ax.add_patch(plt.Rectangle((0.02, y - 0.12), 0.6, 0.11, fill=False, lw=1, ec=INK))
    ax.text(0.32, y - 0.04, t, ha="center", va="center", fontsize=7.5)
    ax.text(0.32, y - 0.09, n, ha="center", va="center", fontsize=8.5, weight="bold")
    if i < len(steps) - 1:
        ax.annotate("", (0.32, y - 0.15), (0.32, y - 0.12), arrowprops=dict(arrowstyle="->", color=INK))
    if i < len(steps) - 1 and excl[i]:
        ax.add_patch(plt.Rectangle((0.68, y - 0.165), 0.31, 0.07, fill=False, lw=0.8, ec="grey"))
        ax.text(0.835, y - 0.13, "Excluded: " + excl[i], ha="center", va="center", fontsize=6.8, color=INK)
        ax.annotate("", (0.68, y - 0.13), (0.32, y - 0.13), arrowprops=dict(arrowstyle="->", color="grey", lw=0.8))
ax.set_xlim(0, 1); ax.set_ylim(-0.1, 1.0); save(fig, "Fig1_flow")

# Fig 2: forest plot
d = pd.read_csv(f"{B}/results/disproportionality.csv")
order = ["main: men, PS", "men 18-39, PS", "men 40-64, PS", "men 65-119, PS", "sens: PS+SS", "sens: no lawyer reports",
         "sens: excl 2014-16", "sens: no duplicate collapse", "active comparator vs PDE5i"]
labels = ["All men (primary analysis)", "  Age 18–39", "  Age 40–64", "  Age ≥65", "Primary or secondary suspect",
          "Lawyer-submitted reports excluded", "First received 2014–2016 excluded", "No exact-duplicate collapse", "Active comparator: PDE5 inhibitors"]
fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2), sharey=True)
for ax, ev, col, title in [(axes[0], "sleep_apnoea", SA, "Sleep apnoea"), (axes[1], "polycythaemia", PC, "Polycythaemia")]:
    sub = d[d.event == ev].set_index("analysis").loc[order]
    y = np.arange(len(order))[::-1]
    ax.errorbar(sub.ROR, y, xerr=[sub.ROR - sub.ROR_lo, sub.ROR_hi - sub.ROR], fmt="s", color=col, ms=4, capsize=2, lw=1)
    for yi, (_, r) in zip(y, sub.iterrows()):
        f = (lambda v: f"{v:.2f}") if r.ROR < 10 else (lambda v: f"{v:.1f}")
        ax.text(1.02, yi, f"{f(r.ROR)} ({f(r.ROR_lo)}–{f(r.ROR_hi)}); n = {int(r.a)}", transform=ax.get_yaxis_transform(), va="center", fontsize=7.5)
    ax.axvline(1, color="grey", ls="--", lw=0.8); ax.set_xscale("log"); ax.set_title(title, loc="left", weight="bold")
    grid = (1, 2, 5, 10) if sub.ROR_hi.max() < 20 else (1, 3, 10, 30, 100, 300)
    ticks = [t for t in grid if t <= sub.ROR_hi.max() * 1.2]
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.minorticks_off()
    ax.set_xlabel("Reporting odds ratio (95% CI), log scale")
axes[0].set_yticks(np.arange(len(order))[::-1]); axes[0].set_yticklabels(labels)
fig.subplots_adjust(wspace=0.9); save(fig, "Fig2_forest")

# Fig 3: cases by FDA receipt year
fig, ax = plt.subplots(figsize=(6.4, 3))
for ev, col, lab in [("sleep_apnoea", SA, "Sleep apnoea"), ("polycythaemia", PC, "Polycythaemia")]:
    pts = ",".join(f"'{p}'" for p in a.EVENTS[ev])
    y = a.c.execute(f"""select yr, count(*) n, count(*) filter (where occp_cod='LW') lw from base where sex='M' and primaryid in (select * from t_ps)
        and primaryid in (select primaryid from reac where pt in ({pts})) and yr>=2012 group by 1 order by 1""").df()
    y = y.set_index("yr").reindex(range(2012, 2027), fill_value=0)
    full = y.loc[:2025]
    ax.plot(full.index, full.n, "-o", color=col, ms=3, lw=1.2, label=lab)
    ax.plot(full.index, full.n - full.lw, ":", color=col, lw=1, label=f"{lab}, lawyer reports removed")
    ax.plot([2026], [y.n[2026]], "o", mfc="white", mec=col, ms=4)
ax.axvspan(2013.5, 2016.5, color="grey", alpha=0.12, lw=0, label="2014–2016 (peak lawyer reporting)")
ax.set_xlabel("Year of initial FDA receipt (open markers: 2026, Q1–Q2 only)"); ax.set_ylabel("Cases (men, testosterone PS)"); ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
save(fig, "Fig3_by_year")

# Fig 4: time to onset
fig, ax = plt.subplots(figsize=(5, 3))
for ev, col, lab in [("sleep_apnoea", SA, "Sleep apnoea"), ("polycythaemia", PC, "Polycythaemia")]:
    t = pd.read_csv(f"{B}/results/tto_{ev}.csv").tto_days.sort_values().reset_index(drop=True)
    m = (t / 30.44).to_numpy(); cum = np.arange(1, len(t) + 1) / len(t); keep = m <= 60
    xs = np.r_[0, m[keep], 60]; ys = np.r_[0, cum[keep], cum[keep][-1] if keep.any() else 0]
    ax.step(xs, ys, where="post", color=col, label=f"{lab} (n = {len(t)}; median {t.median():.0f} days)")
ax.set_xlabel("Months from testosterone start to event"); ax.set_ylabel("Cumulative proportion of cases"); ax.set_xlim(0, 62); ax.set_xticks(range(0, 61, 12))
ax.legend(frameon=False, fontsize=7, loc="lower right"); save(fig, "Fig4_time_to_onset")
print("figures:", sorted(os.listdir(F)))
